#!/usr/bin/env python3
"""Build the two truth read sets of the signal calibration (issue #28).

hg38:    pseudo-reads cut from GRCh38 with samtools faidx at pinned,
         seed-drawn regions (chr1-22, X; 2-8 kb, either strand, no N),
         uppercased. With --regions an existing regions.tsv is re-extracted
         instead of drawing new ones.
uniform: the same number and lengths of reads, cut in order from the
         bases of `phantom-science generate --platform ont --stdout`
         (uniform random composition), so the two sets differ only in
         composition.

copies:  N copies of each record of a FASTA (the RNA004 control
         constructs), interleaved, U written as T.

Read names are name-based UUIDs (uuid5 of the region or chunk label), as
BLOW5/pod5/Dorado require. The reads carry no simulated errors: they are the
truth that signal and squigulator synthesise from.

Usage:
  make_inputs.py hg38 --fasta hg38.fa --out hg38.fa --regions-out regions.tsv \
      [--n 2000] [--seed 28] [--regions regions.tsv]
  make_inputs.py uniform --like hg38.fa --fastq generate.fastq[.gz] --out uniform.fa
  make_inputs.py copies --fasta sampled_context_strands.fa --copies 100 --out rna.fa

Standard library only; run with python3 -I.
"""
import argparse
import gzip
import hashlib
import random
import subprocess
import sys
import uuid

NAMESPACE = uuid.UUID("5f1d0a52-28c0-4d8e-9b1a-7a3f00000028")
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def write_fasta(path, recs):
    h = hashlib.sha256()
    with open(path, "w") as f:
        for name, seq in recs:
            rec = f">{name}\n{seq}\n"
            f.write(rec)
            h.update(rec.encode())
    return h.hexdigest()


def faidx(fasta, regions):
    """Extract regions [(chrom, start1, end1, strand)] in order."""
    out = []
    for chrom, s, e, strand in regions:
        cmd = ["samtools", "faidx", fasta, f"{chrom}:{s}-{e}"]
        if strand == "-":
            cmd.insert(2, "-i")
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        seq = "".join(l for l in res.stdout.splitlines()[1:])
        out.append(seq.upper())
    return out


def hg38(a):
    if a.regions:
        regions = []
        with open(a.regions) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                name, chrom, s, e, strand = line.split()[:5]
                regions.append((name, chrom, int(s), int(e), strand))
        seqs = faidx(a.fasta, [r[1:] for r in regions])
        bad = [r[0] for r, s in zip(regions, seqs) if len(s) != r[3] - r[2] + 1 or "N" in s]
        if bad:
            sys.exit(f"error: {len(bad)} pinned regions changed (first {bad[0]}): wrong reference?")
    else:
        lengths = {}
        with open(a.fasta + ".fai") as f:
            for line in f:
                c, n = line.split()[:2]
                lengths[c] = int(n)
        chroms = [c for c in CHROMS if c in lengths]
        weights = [lengths[c] for c in chroms]
        rng = random.Random(a.seed)
        regions, seqs = [], []
        while len(regions) < a.n:
            cand = []
            for _ in range(a.n - len(regions)):
                c = rng.choices(chroms, weights)[0]
                n = rng.randint(2000, 8000)
                s = rng.randint(1, lengths[c] - n)
                cand.append((c, s, s + n - 1, "-" if rng.random() < 0.5 else "+"))
            for r, seq in zip(cand, faidx(a.fasta, cand)):
                if "N" in seq or len(seq) != r[2] - r[1] + 1:
                    continue
                label = f"{r[0]}:{r[1]}-{r[2]}({r[3]})"
                regions.append((str(uuid.uuid5(NAMESPACE, label)),) + r)
                seqs.append(seq)
        with open(a.regions_out, "w") as f:
            f.write(f"# hg38 pseudo-reads for #28: seed {a.seed}, n {a.n}, 1-based inclusive\n")
            f.write("#name\tchrom\tstart\tend\tstrand\n")
            for r in regions:
                f.write("\t".join(map(str, r)) + "\n")
    digest = write_fasta(a.out, [(r[0], s) for r, s in zip(regions, seqs)])
    gc = sum(s.count("G") + s.count("C") for s in seqs) / sum(map(len, seqs))
    print(f"hg38\t{len(seqs)} reads\t{sum(map(len, seqs))} bases\tGC {gc:.4f}\tsha256 {digest}")


def uniform(a):
    lengths = [len(s) for _, s in read_fasta(a.like)]
    pool, need = [], sum(lengths)
    have = 0
    with open(a.fastq, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    with (gzip.open(a.fastq, "rt") if gz else open(a.fastq)) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                s = line.strip().upper()
                pool.append(s)
                have += len(s)
                if have >= need:
                    break
    bases = "".join(pool)
    if len(bases) < need:
        sys.exit(f"error: generate output has {len(bases)} bases, need {need}")
    recs, at = [], 0
    for i, n in enumerate(lengths):
        recs.append((str(uuid.uuid5(NAMESPACE, f"uniform:{i}")), bases[at:at + n]))
        at += n
    digest = write_fasta(a.out, recs)
    seq = bases[:need]
    gc = (seq.count("G") + seq.count("C")) / need
    print(f"uniform\t{len(recs)} reads\t{need} bases\tGC {gc:.4f}\tsha256 {digest}")


def copies(a):
    recs = []
    for name, seq in read_fasta(a.fasta):
        for i in range(a.copies):
            recs.append((str(uuid.uuid5(NAMESPACE, f"{name}:{i}")), seq.upper().replace("U", "T")))
    # Interleave constructs so any prefix of the file covers all of them.
    n = len(recs) // a.copies
    recs = [recs[j * a.copies + i] for i in range(a.copies) for j in range(n)]
    digest = write_fasta(a.out, recs)
    print(f"copies\t{len(recs)} reads\t{n} constructs x {a.copies}\tsha256 {digest}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("copies")
    c.add_argument("--fasta", required=True)
    c.add_argument("--copies", type=int, default=100)
    c.add_argument("--out", required=True)
    h = sub.add_parser("hg38")
    h.add_argument("--fasta", required=True)
    h.add_argument("--out", required=True)
    h.add_argument("--regions-out", default="regions.tsv")
    h.add_argument("--regions")
    h.add_argument("--n", type=int, default=2000)
    h.add_argument("--seed", type=int, default=28)
    u = sub.add_parser("uniform")
    u.add_argument("--like", required=True)
    u.add_argument("--fastq", required=True)
    u.add_argument("--out", required=True)
    a = p.parse_args()
    {"hg38": hg38, "uniform": uniform, "copies": copies}[a.cmd](a)


if __name__ == "__main__":
    main()
