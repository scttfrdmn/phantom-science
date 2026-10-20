#!/usr/bin/env python3
"""Does per-sample noise depend on the k-mer? (issue #28)

Segments real (or synthetic) signal into bases with Dorado's moves (see
measure.py), assigns each base to the called k-mer at --offset, and
estimates each k-mer's per-sample noise from the within-base first
differences (mean |Δx| / (2/√π) / √2 over bases with at least 4 samples,
so a level step at a misplaced move boundary inflates it equally for all
k-mers). Reports the correlation of that per-k-mer noise with the level and
with the level_stdv column of --stdv-table (e.g. squigulator's built-in
table recovered with squig_table.py), over k-mers with at least --min-count
bases, and a split-half correlation of the noise itself (its reliability).

Usage:
  kmer_noise.py --sam calls.sam --pod5 reads.pod5 --stdv-table T.txt \
      [--offset 6] [--rna] [--max-reads N]

Run with the pod5/numpy venv's python, -I.
"""
import argparse
import os
import sys

import numpy as np
import pod5

sys.dont_write_bytecode = True  # no __pycache__ next to the scripts
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measure import CODE, kmer_index, parse_sam  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sam", required=True)
    p.add_argument("--pod5", required=True)
    p.add_argument("--stdv-table", required=True)
    p.add_argument("--offset", type=int, default=6)
    p.add_argument("--rna", action="store_true")
    p.add_argument("--max-reads", type=int, default=100000)
    p.add_argument("--min-count", type=int, default=20)
    a = p.parse_args()

    rows = [l.split() for l in open(a.stdv_table) if l.strip() and not l.startswith("#")]
    rows = [r for r in rows if r[0].upper() != "KMER"]
    k = len(rows[0][0])
    level = np.zeros(4 ** k)
    stdv = np.zeros(4 ** k)
    for r in rows:
        idx = 0
        for b in r[0].encode():
            idx = idx * 4 + CODE[b]
        level[idx], stdv[idx] = float(r[1]), float(r[2])

    reads = {}
    for rid, seq, ts, stride, mv, split in parse_sam(a.sam):
        starts = np.flatnonzero(mv == 1)
        if split or len(starts) < 20 or len(seq) != len(starts):
            continue
        reads[rid] = (ts + starts * stride, ts + stride * len(mv), seq)
        if len(reads) >= a.max_reads:
            break
    nk = 4 ** k
    s = [np.zeros(nk), np.zeros(nk)]
    c = [np.zeros(nk), np.zeros(nk)]
    with pod5.Reader(a.pod5) as r:
        for j, rec in enumerate(r.reads(selection=list(reads), missing_ok=True)):
            bounds, end, seq = reads[str(rec.read_id)]
            x = rec.signal_pa.astype(np.float64)
            if bounds[-1] >= len(x):
                continue
            edges = np.append(bounds, min(end, len(x)))
            lens = np.diff(edges)
            if (lens <= 0).any():
                continue
            xs = x[edges[0]:edges[-1]]
            rel = edges[:-1] - edges[0]
            seg = np.repeat(np.arange(len(lens)), lens)
            ad = np.abs(np.diff(xs))
            same = seg[:-1] == seg[1:]
            sums = np.bincount(seg[:-1][same], weights=ad[same], minlength=len(lens))
            cnts = np.bincount(seg[:-1][same], minlength=len(lens))
            n = len(seq)
            pos = np.arange(n)[::-1] if a.rna else np.arange(n)
            kidx = kmer_index(seq.upper(), k)
            st = pos - a.offset
            ok = (st >= 0) & (st < n) & (lens >= 4)
            ok[ok] &= kidx[st[ok]] >= 0
            half = j % 2
            ki = kidx[st[ok]]
            np.add.at(s[half], ki, sums[ok] / cnts[ok])
            np.add.at(c[half], ki, 1)
            del rel
    conv = np.sqrt(np.pi) / 2 / np.sqrt(2)
    tot_s, tot_c = s[0] + s[1], c[0] + c[1]
    ok = tot_c >= a.min_count
    noise = np.where(ok, tot_s / np.maximum(tot_c, 1) * conv, np.nan)
    both = (c[0] >= a.min_count / 2) & (c[1] >= a.min_count / 2)
    h0 = s[0][both] / c[0][both]
    h1 = s[1][both] / c[1][both]
    print(f"reads {len(reads)}; k-mers with >= {a.min_count} bases: {ok.sum()}")
    print(f"per-k-mer noise (pA): median {np.nanmedian(noise):.2f}, 5-95% "
          f"{np.nanpercentile(noise, 5):.2f}-{np.nanpercentile(noise, 95):.2f}")
    print(f"split-half reliability r {np.corrcoef(h0, h1)[0, 1]:.3f} over {both.sum()} k-mers")
    print(f"corr(noise, level) {np.corrcoef(noise[ok], level[ok])[0, 1]:.3f}")
    print(f"corr(noise, table level_stdv) {np.corrcoef(noise[ok], stdv[ok])[0, 1]:.3f}  "
          f"(table stdv median {np.median(stdv):.2f} pA)")
    b1, b0 = np.polyfit(stdv[ok], noise[ok], 1)
    print(f"noise = {b0:.2f} + {b1:.3f} x table stdv")


if __name__ == "__main__":
    sys.exit(main())
