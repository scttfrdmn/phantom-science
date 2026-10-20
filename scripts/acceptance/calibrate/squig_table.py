#!/usr/bin/env python3
"""Recover squigulator's built-in k-mer level table (issue #28).

squigulator embeds its pore models in the binary. This script writes a de
Bruijn sequence B(4, k), in which every k-mer occurs once, for squigulator
to simulate, then reads the levels back:

  debruijn.py-style:  squig_table.py debruijn --k 9 --out db.fa
  squigulator -x dna-r10-prom --full-contigs --ideal -o ideal.blow5 db.fa
  squigulator -x dna-r10-prom --full-contigs --ideal-time -o amp.blow5 db.fa
  squig_table.py levels --ideal ideal.blow5 --amp amp.blow5 --fasta db.fa \
      --k 9 --out squig_pa.txt [--compare ont_table --scale S --shift H]

With --ideal (no amplitude or time noise) each k-mer is a constant run of
samples; runs are split at level changes and assigned to the k-mers in
order. With --ideal-time (amplitude noise only, same dwells) the SD of each
run is that k-mer's per-sample stdv. The output is a pA table with a
level_stdv column (signal reads it as pA). --compare reports the
correlation and RMS difference with another table.

Run with the pyslow5/numpy venv's python, -I.
"""
import argparse
import sys

import numpy as np


def debruijn(k, n=4):
    a = [0] * (k * n)
    seq = []

    def db(t, p):
        if t > k:
            if k % p == 0:
                seq.extend(a[1:p + 1])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, n):
                a[t] = j
                db(t + 1, t)

    sys.setrecursionlimit(10000)
    db(1, 1)
    s = "".join("ACGT"[i] for i in seq)
    return s + s[:k - 1]


def read_one(path):
    import pyslow5
    f = pyslow5.Open(path, "r")
    for rec in f.seq_reads(pA=True):
        return np.asarray(rec["signal"], dtype=np.float64)
    raise SystemExit(f"no reads in {path}")


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("debruijn")
    d.add_argument("--k", type=int, default=9)
    d.add_argument("--out", required=True)
    lv = sub.add_parser("levels")
    lv.add_argument("--ideal", required=True)
    lv.add_argument("--amp")
    lv.add_argument("--fasta", required=True)
    lv.add_argument("--k", type=int, default=9)
    lv.add_argument("--out", required=True)
    lv.add_argument("--compare")
    lv.add_argument("--scale", type=float, default=1.0)
    lv.add_argument("--shift", type=float, default=0.0)
    lv.add_argument("--reverse", action="store_true",
                    help="the signal walks the sequence 3'->5' (direct RNA)")
    a = p.parse_args()

    if a.cmd == "debruijn":
        s = debruijn(a.k)
        with open(a.out, "w") as f:
            f.write(f">debruijn_{a.k}\n{s}\n")
        print(f"wrote {a.out}: {len(s)} bases")
        return

    seq = "".join(l.strip() for l in open(a.fasta) if not l.startswith(">"))
    nk = len(seq) - a.k + 1
    x = read_one(a.ideal)
    change = np.flatnonzero(np.diff(x) != 0) + 1
    starts = np.concatenate([[0], change])
    lens = np.diff(np.append(starts, len(x)))
    print(f"ideal signal: {len(x)} samples, {len(starts)} runs, {nk} k-mers")
    # Equal consecutive levels merge two runs; a run twice the usual length
    # is split in half (rare; reported).
    mean_len = np.median(lens)
    runs = []
    for s0, ln in zip(starts, lens):
        parts = max(1, int(round(ln / mean_len)))
        for j in range(parts):
            runs.append((s0 + j * ln // parts, ln // parts))
    if len(runs) != nk:
        raise SystemExit(f"error: {len(runs)} runs after splitting, want {nk}; "
                         "is the profile's prefix off (--prefix=no) and the read forward?")
    amp = read_one(a.amp) if a.amp else None
    if amp is not None and len(amp) != len(x):
        raise SystemExit("error: --amp signal length differs; dwells must be identical (--ideal-time)")
    code = {c: i for i, c in enumerate("ACGT")}
    rows = {}
    if a.reverse:
        runs = runs[::-1]
    for i, (s0, ln) in enumerate(runs):
        km = seq[i:i + a.k]
        idx = 0
        for c in km:
            idx = idx * 4 + code[c]
        mean = x[s0:s0 + ln].mean()
        sd = amp[s0:s0 + ln].std(ddof=1) if amp is not None and ln > 2 else float("nan")
        rows[idx] = (km, mean, sd)
    if len(rows) != 4 ** a.k:
        raise SystemExit(f"error: {len(rows)} distinct k-mers, want {4 ** a.k}")
    sds = np.array([r[2] for r in rows.values()])
    fill = np.nanmedian(sds) if amp is not None else 0.08 * 23.506
    with open(a.out, "w") as f:
        f.write(f"#k\t{a.k}\nkmer\tlevel_mean\tlevel_stdv\n")
        for idx in sorted(rows):
            km, mean, sd = rows[idx]
            f.write(f"{km}\t{mean:.4f}\t{(sd if np.isfinite(sd) else fill):.4f}\n")
    lvl = np.array([rows[i][1] for i in sorted(rows)])
    print(f"wrote {a.out}: level mean {lvl.mean():.2f} SD {lvl.std():.2f} pA; "
          f"stdv median {np.nanmedian(sds):.2f} pA (5-95%: {np.nanpercentile(sds, 5):.2f}-{np.nanpercentile(sds, 95):.2f})")
    if a.compare:
        other = np.zeros(4 ** a.k)
        for l in open(a.compare):
            c = l.split()
            if not c or c[0].startswith("#") or c[0].lower() == "kmer":
                continue
            idx = 0
            for ch in c[0]:
                idx = idx * 4 + code[ch]
            other[idx] = float(c[1]) * a.scale + a.shift
        r = np.corrcoef(lvl, other)[0, 1]
        b1, b0 = np.polyfit(other, lvl, 1)
        print(f"vs {a.compare}: r {r:.5f}; squigulator = {b0:.3f} + {b1:.4f} x other; "
              f"RMS difference {np.sqrt(np.mean((lvl - other) ** 2)):.2f} pA, "
              f"after the linear fit {np.std(lvl - (b0 + b1 * other)):.2f} pA")


if __name__ == "__main__":
    sys.exit(main())
