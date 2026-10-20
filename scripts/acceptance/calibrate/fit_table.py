#!/usr/bin/env python3
"""Compare the k-mer level table with real current, and refit it (issue #28).

Both real and synthetic signal are segmented into bases by Dorado's move
table (see measure.py). For each base, its mean level is assigned to the
called k-mer at --offset (the offset where measure.py's table fit peaks)
and normalised per read by a linear fit against the table, so per-read
shift and scale drop out. Averaging over occurrences gives an empirical
level per k-mer: E_real from real reads, E_syn from signal synthesised
from the table itself.

Dorado's segmentation blurs levels the same way for both, so E_syn is the
table as seen through that measurement and E_real − E_syn estimates the
table's error as seen through it. The script reports how well the table
describes each (correlations, attenuation slope of E_syn on the table, SD
of the shrunken difference, CpG vs non-CpG mean difference) and, with
--out, writes a corrected two-column table in the input table's units:

  T' = T + w·(E_real − E_syn) / slope

where w is the empirical-Bayes shrinkage of each k-mer's difference towards
0 (from its read counts) and slope the attenuation. With --exclude-cpg,
k-mers containing CG keep their table level (native human DNA is CpG
methylated; the synthetic signal models canonical bases).

Usage:
  fit_table.py --real-sam real.sam --real-pod5 real.pod5 \
      --syn-sam syn.sam --syn-pod5 syn.pod5 --table 9mer_levels_v1.txt \
      --scale 23.506 --shift 93.692 [--offset 6] [--exclude-cpg] [--out T.txt] \
      [--stdv-out T3.txt]

It also prints the per-sample noise of real and synthetic signal by table
level; --stdv-out writes the table in pA with a level_stdv column scaled by
the real/synthetic noise ratio of each level bin.

Caveat (results/dna.tsv, dna_hg38_fittedtable): the refit table lowered
identity from 91.75% to 77.54%. Dorado's move boundaries are not level
boundaries, so these per-k-mer means are not a substitute for resquiggling.

Run with the pod5/numpy venv's python, -I.
"""
import argparse
import os
import sys

import numpy as np
import pod5

sys.dont_write_bytecode = True  # no __pycache__ next to the scripts
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measure import kmer_index, load_table, parse_sam, robust_sd  # noqa: E402


def accumulate(sam, pod5_path, k, table, offset, max_reads, rna):
    nk = len(table)
    s = np.zeros(nk)
    ss = np.zeros(nk)
    c = np.zeros(nk)
    reads = {}
    for rid, seq, ts, stride, mv, split in parse_sam(sam):
        starts = np.flatnonzero(mv == 1)
        if split or len(starts) < 50 or len(seq) != len(starts):
            continue
        reads[rid] = (ts + starts * stride, ts + stride * len(mv), seq)
        if len(reads) >= max_reads:
            break
    resid_ss = resid_n = 0.0
    noise_lv, noise_dx = [], []
    with pod5.Reader(pod5_path) as r:
        for rec in r.reads(selection=list(reads), missing_ok=True):
            bounds, end, seq = reads[str(rec.read_id)]
            x = rec.signal_pa.astype(np.float64)
            if bounds[-1] >= len(x):
                continue
            edges = np.append(bounds, min(end, len(x)))
            lens = np.diff(edges)
            if (lens <= 0).any():
                continue
            m = np.add.reduceat(x[edges[0]:edges[-1]], edges[:-1] - edges[0]) / lens
            n = len(seq)
            pos = np.arange(n)[::-1] if rna else np.arange(n)
            kidx = kmer_index(seq.upper(), k)
            st = pos - offset
            ok = (st >= 0) & (st < n)
            ok[ok] &= kidx[st[ok]] >= 0
            if ok.sum() < 50:
                continue
            ki = kidx[st[ok]]
            lv, mm = table[ki], m[ok]
            if lv.std() == 0 or mm.std() == 0:
                continue
            b1, b0 = np.polyfit(lv, mm, 1)
            if not b1 > 0.2:
                continue
            norm = (mm - b0) / b1
            res = norm - lv
            sd = res.std()
            keep = np.abs(res) < 4 * sd
            resid_ss += float(np.sum(res[keep] ** 2))
            resid_n += int(keep.sum())
            np.add.at(s, ki[keep], norm[keep])
            np.add.at(ss, ki[keep], norm[keep] ** 2)
            np.add.at(c, ki[keep], 1)
            if len(noise_lv) < NOISE_READS:
                # Within-base |Δx| tagged with the base's table level.
                seg = np.repeat(np.arange(len(lens)), lens)
                xs = x[edges[0]:edges[-1]]
                same = seg[:-1] == seg[1:]
                lvl = np.full(len(lens), np.nan)
                lvl[ok] = lv
                noise_lv.append(lvl[seg[:-1][same]])
                noise_dx.append(np.abs(np.diff(xs))[same])
    lv_all = np.concatenate(noise_lv)
    dx_all = np.concatenate(noise_dx)
    good = ~np.isnan(lv_all)
    return s, ss, c, np.sqrt(resid_ss / max(resid_n, 1)), len(reads), (lv_all[good], dx_all[good])


NOISE_READS = 300
BINS = np.arange(40, 151, 10)


def binned_noise(lv, dx):
    """Robust per-sample noise SD (measure.robust_sd) per level bin."""
    out = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        sel = (lv >= lo) & (lv < hi)
        out.append(robust_sd(dx[sel]) if sel.sum() > 1000 else np.nan)
    return np.array(out)


def binned_all(noise):
    return robust_sd(noise[1])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--real-sam", required=True)
    p.add_argument("--real-pod5", required=True)
    p.add_argument("--syn-sam", required=True)
    p.add_argument("--syn-pod5", required=True)
    p.add_argument("--table", required=True)
    p.add_argument("--scale", type=float, default=1.0)
    p.add_argument("--shift", type=float, default=0.0)
    p.add_argument("--offset", type=int, default=6)
    p.add_argument("--max-reads", type=int, default=100000)
    p.add_argument("--min-count", type=int, default=3)
    p.add_argument("--exclude-cpg", action="store_true")
    p.add_argument("--rna", action="store_true")
    p.add_argument("--out")
    p.add_argument("--stdv-out", help="write a pA table with a per-k-mer level_stdv column")
    p.add_argument("--syn-stdv", type=float, default=0.08 * 23.506,
                   help="per-sample stdv (pA) the synthetic signal was made with")
    a = p.parse_args()

    k, T = load_table(a.table, a.scale, a.shift)
    sr, ssr, cr, rr, nr, noise_r = accumulate(a.real_sam, a.real_pod5, k, T, a.offset, a.max_reads, a.rna)
    sy, ssy, cy, ry, ny, noise_y = accumulate(a.syn_sam, a.syn_pod5, k, T, a.offset, a.max_reads, a.rna)
    nb_r, nb_y = binned_noise(*noise_r), binned_noise(*noise_y)
    print("per-sample noise SD by table level (pA; robust, within Dorado bases):")
    print("  level   " + " ".join(f"{lo:>6d}" for lo in BINS[:-1]))
    print("  real    " + " ".join(f"{v:6.2f}" for v in nb_r))
    print("  synth   " + " ".join(f"{v:6.2f}" for v in nb_y))
    print(f"  overall real {binned_all(noise_r):.2f} synthetic {binned_all(noise_y):.2f}")
    both = (cr >= a.min_count) & (cy >= a.min_count)
    Er = np.where(cr > 0, sr / np.maximum(cr, 1), np.nan)
    Ey = np.where(cy > 0, sy / np.maximum(cy, 1), np.nan)
    # Per-observation residual variance around each k-mer mean, pooled.
    vr = np.nansum(ssr - cr * np.nan_to_num(Er) ** 2) / max(cr.sum() - (cr > 0).sum(), 1)
    vy = np.nansum(ssy - cy * np.nan_to_num(Ey) ** 2) / max(cy.sum() - (cy > 0).sum(), 1)
    slope = np.polyfit(T[both], Ey[both], 1)[0]
    d = np.zeros_like(T)
    d[both] = Er[both] - Ey[both]
    noise = np.zeros_like(T)
    noise[both] = vr / cr[both] + vy / cy[both]
    tau2 = max(np.var(d[both]) - np.mean(noise[both]), 0.0)
    w = np.zeros_like(T)
    w[both] = tau2 / (tau2 + noise[both])
    names = [None] * len(T)
    cpg = np.zeros(len(T), dtype=bool)
    for i in range(len(T)):
        s = "".join("ACGT"[(i >> (2 * (k - 1 - j))) & 3] for j in range(k))
        names[i] = s
        cpg[i] = "CG" in s
    print(f"reads real {nr} synthetic {ny}; k-mers with >= {a.min_count} obs in both: {both.sum()}")
    print(f"per-base residual SD around the table (pA, per-read fit): real {rr:.2f} synthetic {ry:.2f}")
    print(f"per-observation SD around the k-mer mean (pA): real {np.sqrt(vr):.2f} synthetic {np.sqrt(vy):.2f}")
    print(f"corr(E_real, table) {np.corrcoef(Er[both], T[both])[0, 1]:.4f}  "
          f"corr(E_syn, table) {np.corrcoef(Ey[both], T[both])[0, 1]:.4f}  "
          f"corr(E_real, E_syn) {np.corrcoef(Er[both], Ey[both])[0, 1]:.4f}")
    print(f"attenuation slope of E_syn on table {slope:.3f}")
    print(f"SD of E_real - E_syn {np.std(d[both]):.2f} pA; noise part {np.sqrt(np.mean(noise[both])):.2f}; "
          f"signal part (tau) {np.sqrt(tau2):.2f} pA -> table error ~ {np.sqrt(tau2) / slope:.2f} pA")
    for lab, sel in (("CpG", cpg & both), ("non-CpG", ~cpg & both)):
        print(f"  {lab:8s} k-mers {sel.sum():6d}  mean E_real - E_syn {d[sel].mean():+.2f} pA  SD {d[sel].std():.2f}")
    if a.out:
        corr = w * d / slope
        if a.exclude_cpg:
            corr[cpg] = 0
        new = (T + corr - a.shift) / a.scale
        with open(a.out, "w") as f:
            for i in range(len(T)):
                f.write(f"{names[i]}\t{new[i]:.6f}\n")
        print(f"wrote {a.out}: mean |correction| {np.abs(corr).mean():.2f} pA, "
              f"{(corr != 0).sum()} k-mers changed")
    if a.stdv_out:
        # Per-k-mer stdv in pA: the synthetic run's stdv (--syn-stdv, the
        # constant DefaultModelStdv x scale it was made with) times the
        # real/synthetic noise ratio of the k-mer's level bin.
        ratio = nb_r / nb_y
        fill = np.nanmedian(ratio)
        ratio = np.where(np.isnan(ratio), fill, ratio)
        b = np.clip(np.digitize(T, BINS) - 1, 0, len(ratio) - 1)
        sd = a.syn_stdv * ratio[b]
        with open(a.stdv_out, "w") as f:
            f.write("#k\t%d\nkmer\tlevel_mean\tlevel_stdv\n" % k)
            for i in range(len(T)):
                f.write(f"{names[i]}\t{T[i]:.4f}\t{sd[i]:.4f}\n")
        print(f"wrote {a.stdv_out}: stdv {sd.min():.2f}-{sd.max():.2f} pA (ratio by bin "
              + " ".join(f"{r:.2f}" for r in ratio) + ")")


if __name__ == "__main__":
    sys.exit(main())
