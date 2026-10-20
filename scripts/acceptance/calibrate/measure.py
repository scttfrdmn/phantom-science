#!/usr/bin/env python3
"""Score one calibration condition (issue #28).

Reads Dorado calls (SAM text from `samtools view`, with move tables), the
minimap2 PAF of those calls and the pod5 the calls were made from, and
prints one key=value per line.

identity (PAF primary alignments, matches / block length, as
dorado_identity.sh):
  n_reads n_aligned median_identity ci_lo ci_hi mean_identity
    ci_* is the percentile bootstrap 95% CI of the median over reads
    (--boot resamples, numpy seed --seed). --min-mapq filters alignments
    (20 for real reads against GRCh38).
  aligned_frac   aligned reads / basecalled reads
  median_tident  median over basecalled reads of matches / target length
                 (0 if unaligned): penalises truncated and missing calls;
                 meaningful only when the targets are the read truths
                 (short RNA constructs), not a genome

dwell (Dorado moves, mv:B:c,stride,...; the last base of a read is open and
skipped):
  bases mean_dwell dwell_cv
  dwell_short    fraction of bases at the minimum dwell (one stride)

signal (pA from pod5, over the moved samples after ts:i, segmented into
bases by the moves; split reads, pi:Z, are skipped; at most --signal-reads
reads):
  noise_sd       pooled within-base SD, pA (includes level changes Dorado's
                 move boundaries do not line up with)
  noise_mad      segmentation-free per-sample noise, robust to level
                 steps: robust_sd of the |Δx| (pA)
  diff_ac1       lag-1 autocorrelation of Δx: −0.5 for white noise on a
                 step function, larger when levels change smoothly
  seg_mean_sd    SD of per-base mean levels, pA
  step_abs       mean |difference| of consecutive base means, pA
  resid_ac1      lag-1 autocorrelation of within-base residuals
  carry1 carry2  for steps |Δmean| > --step-pa between consecutive bases,
                 mean of (x_i − m) / (m_prev − m) for the first and second
                 sample of the new base
  signal_mean signal_sd

level table (with --table): per-base mean levels against the table level of
the called k-mer, at the k-mer offset that fits best (searched over
−M..k+M−1 bases, M = --max-offset):
  table_r        mean per-read Pearson r
  table_resid_sd pooled residual SD of a per-read linear fit, pA
  table_offset   best offset (bases)
  For synthetic signal this measures how well Dorado's moves and calls
  recover the table it was made from; for real reads, how well the table
  describes real current.

Run with the pod5/numpy venv's python, -I.
"""
import argparse
import sys

import numpy as np
import pod5

CODE = np.full(256, -1, dtype=np.int64)
for i, b in enumerate(b"ACGT"):
    CODE[b] = i
CODE[ord("U")] = 3


def robust_sd(absdx):
    """Per-sample Gaussian SD from first differences |Δx|, ignoring the top
    5% (level steps): for Z ~ N(0,1), E[|Z| ; |Z| < 1.96] / 0.95 = 0.7169,
    and Δx has SD σ·√2. A trimmed mean rather than the median, which the
    ADC step (0.14-0.19 pA) quantises."""
    if len(absdx) == 0:
        return float("nan")
    cut = np.percentile(absdx, 95)
    return float(absdx[absdx <= cut].mean() / 0.7169 / np.sqrt(2))


def parse_paf(path, min_mapq):
    best = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            tags = {t[:5]: t[5:] for t in c[12:]}
            if tags.get("tp:A:") != "P" or c[0] in best:
                continue
            if int(c[11]) < min_mapq:
                continue
            best[c[0]] = (int(c[9]) / int(c[10]), int(c[9]), int(c[6]))
    return best


def parse_sam(path):
    """Yield (read_id, seq, ts, stride, moves, split) per primary record."""
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.rstrip("\n").split("\t")
            if int(c[1]) & 0x900:
                continue
            mv, ts, split = None, 0, False
            for t in c[11:]:
                if t.startswith("mv:B:c,"):
                    mv = np.array(t[7:].split(","), dtype=np.int64)
                elif t.startswith("ts:i:"):
                    ts = int(t[5:])
                elif t.startswith("pi:Z:"):
                    split = True
            if mv is None or len(mv) < 3:
                continue
            yield c[0], c[9], ts, int(mv[0]), mv[1:], split


def load_table(path, scale, shift):
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
    rows = [r for r in rows if r[0].upper() != "KMER"]
    k = len(rows[0][0])
    lv = np.zeros(4 ** k)
    for r in rows:
        idx = 0
        for b in r[0].encode():
            idx = idx * 4 + CODE[b]
        lv[idx] = float(r[1]) * scale + shift
    return k, lv


def kmer_index(seq, k):
    """Index of the k-mer starting at each position (−1 past the end)."""
    b = CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(b) - k + 1
    out = np.full(len(b), -1, dtype=np.int64)
    if n <= 0 or (b < 0).any():
        return out
    idx = np.zeros(n, dtype=np.int64)
    for t in range(k):
        idx = idx * 4 + b[t:t + n]
    out[:n] = idx
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sam", required=True)
    p.add_argument("--paf", required=True)
    p.add_argument("--pod5", required=True)
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument("--boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--signal-reads", type=int, default=400)
    p.add_argument("--step-pa", type=float, default=8.0)
    p.add_argument("--table", help="k-mer level table for the level-table fit")
    p.add_argument("--table-scale", type=float, default=1.0)
    p.add_argument("--table-shift", type=float, default=0.0)
    p.add_argument("--rna", action="store_true", help="signal runs 3'->5' of the called sequence")
    p.add_argument("--max-offset", type=int, default=12)
    a = p.parse_args()

    # Identity.
    ident = parse_paf(a.paf, a.min_mapq)
    ids = np.array(sorted(v[0] for v in ident.values()))
    out = {"n_aligned": len(ids)}
    if len(ids):
        rng = np.random.default_rng(a.seed)
        boot = np.median(rng.choice(ids, size=(a.boot, len(ids)), replace=True), axis=1)
        out.update(median_identity=np.median(ids), ci_lo=np.percentile(boot, 2.5),
                   ci_hi=np.percentile(boot, 97.5), mean_identity=ids.mean())

    # Dwell, and the reads kept for signal statistics.
    dw_all, segs, n_reads = [], {}, 0
    for rid, seq, ts, stride, mv, split in parse_sam(a.sam):
        n_reads += 1
        starts = np.flatnonzero(mv == 1)
        if len(starts) < 3:
            continue
        dw_all.append(np.diff(starts) * stride)
        if not split and len(segs) < a.signal_reads and len(seq) == len(starts):
            segs[rid] = (ts + starts * stride, stride * len(mv) + ts, seq)
    out["n_reads"] = n_reads
    if n_reads:
        tid = [v[1] / v[2] for v in ident.values()] + [0.0] * max(n_reads - len(ident), 0)
        out.update(aligned_frac=len(ident) / n_reads, median_tident=float(np.median(tid)))
    if dw_all:
        dw = np.concatenate(dw_all).astype(float)
        out.update(bases=len(dw), mean_dwell=dw.mean(), dwell_cv=dw.std() / dw.mean(),
                   dwell_short=float(np.mean(dw <= dw.min())))

    table = load_table(a.table, a.table_scale, a.table_shift) if a.table else None
    k = table[0] if table else 0
    offsets = list(range(-a.max_offset, k + a.max_offset))
    fit_r = {o: [] for o in offsets}
    fit_ss = {o: [0.0, 0.0] for o in offsets}

    ss_within = n_within = 0.0
    ac_num = ac_den = 0.0
    dac_num = dac_den = 0.0
    means, steps, carry1, carry2, absdiff = [], [], [], [], []
    sig_sum = sig_ss = sig_n = 0.0
    if segs:
        with pod5.Reader(a.pod5) as r:
            for rec in r.reads(selection=list(segs), missing_ok=True):
                bounds, end, seq = segs[str(rec.read_id)]
                x = rec.signal_pa.astype(np.float64)
                if bounds[-1] >= len(x):
                    continue
                edges = np.append(bounds, min(end, len(x)))
                x = x[edges[0]:edges[-1]]
                rel = edges[:-1] - edges[0]
                lens = np.diff(edges)
                if (lens <= 0).any():
                    continue
                m = np.add.reduceat(x, rel) / lens
                seg_id = np.repeat(np.arange(len(lens)), lens)
                resid = x - m[seg_id]
                ss_within += float(np.sum(resid ** 2))
                n_within += float(len(x) - len(lens))
                same = seg_id[:-1] == seg_id[1:]
                ac_num += float(np.sum(resid[:-1][same] * resid[1:][same]))
                ac_den += float(np.sum(resid[:-1][same] ** 2))
                dx = np.diff(x)
                absdiff.append(np.abs(dx))
                dxc = dx - dx.mean()
                dac_num += float(np.sum(dxc[:-1] * dxc[1:]))
                dac_den += float(np.sum(dxc ** 2))
                means.append(m)
                steps.append(np.abs(np.diff(m)))
                dm = m[:-1] - m[1:]  # previous minus current
                big = np.flatnonzero(np.abs(dm) > a.step_pa) + 1
                carry1.append((x[rel[big]] - m[big]) / dm[big - 1])
                two = big[lens[big] >= 2]
                carry2.append((x[rel[two] + 1] - m[two]) / dm[two - 1])
                sig_sum += float(x.sum())
                sig_ss += float(np.sum(x ** 2))
                sig_n += len(x)
                if table:
                    # Base i of the signal is called-sequence position pos[i].
                    n = len(seq)
                    pos = np.arange(n)[::-1] if a.rna else np.arange(n)
                    kidx = kmer_index(seq.upper(), k)
                    for o in offsets:
                        st = pos - o
                        ok = (st >= 0) & (st < n)
                        ok[ok] &= kidx[st[ok]] >= 0
                        if ok.sum() < 50:
                            continue
                        lv = table[1][kidx[st[ok]]]
                        mm = m[ok]
                        if lv.std() == 0 or mm.std() == 0:
                            continue
                        rr = np.corrcoef(lv, mm)[0, 1]
                        b1, b0 = np.polyfit(lv, mm, 1)
                        res = mm - (b0 + b1 * lv)
                        fit_r[o].append(rr)
                        fit_ss[o][0] += float(np.sum(res ** 2))
                        fit_ss[o][1] += len(res) - 2
    if sig_n:
        mu = sig_sum / sig_n
        ad = np.concatenate(absdiff)
        out.update(noise_sd=np.sqrt(ss_within / n_within),
                   noise_mad=robust_sd(ad), diff_ac1=dac_num / dac_den,
                   seg_mean_sd=np.concatenate(means).std(), step_abs=np.concatenate(steps).mean(),
                   resid_ac1=ac_num / ac_den,
                   carry1=np.concatenate(carry1).mean(), carry2=np.concatenate(carry2).mean(),
                   signal_mean=mu, signal_sd=np.sqrt(sig_ss / sig_n - mu * mu),
                   signal_reads=len(means))
    if table and any(fit_r.values()):
        best = max((o for o in offsets if fit_r[o]), key=lambda o: np.mean(fit_r[o]))
        out.update(table_r=float(np.mean(fit_r[best])),
                   table_resid_sd=float(np.sqrt(fit_ss[best][0] / fit_ss[best][1])),
                   table_offset=best,
                   table_r_by_offset=",".join(f"{o}:{np.mean(fit_r[o]):.3f}" for o in offsets if fit_r[o]))
    for key, v in out.items():
        print(f"{key}={v:.6g}" if isinstance(v, float) else f"{key}={v}")


if __name__ == "__main__":
    sys.exit(main())
