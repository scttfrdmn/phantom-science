"""Measure Illumina quality, error and insert-size statistics (#27).

Run with `python3 -I` (no third-party packages). Subcommands:

  fastq OUTDIR R1.fastq R2.fastq
      Quality alphabet, per-cycle quality frequencies, per-cycle transition
      counts P(q[c+1] | q[c]), run lengths, low-quality bases per read and
      N calls, for each mate.

  sam OUTDIR < name-grouped SAM (minimap2 output, with --MD)
      Substitution rate per quality and per cycle (primary alignments with
      MAPQ >= 30; mismatches from MD, indels and soft-clipped bases
      excluded; N bases and reference N excluded), indel rates, the
      mate-overlap error estimate (a mismatch in one mate where the other
      mate agrees with the reference is a sequencing error; the same
      non-reference base in both mates is a variant of the sample or a
      pre-sequencing change) and insert sizes of proper pairs.

All output is tab-separated with a header line.
"""

import itertools
import math
import os
import re
import statistics
import sys

MIN_MAPQ = 30
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def write_tsv(path, header, rows):
    with open(path, "w", encoding="ascii") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def fmt(x, digits=6):
    return f"{x:.{digits}g}"


# ---------------------------------------------------------------- fastq


def fastq_records(path):
    with open(path, "r", encoding="ascii") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip("\n")
            f.readline()
            q = f.readline().rstrip("\n")
            if len(s) != len(q):
                raise SystemExit(f"{path}: sequence and quality lengths differ at {h.strip()}")
            yield s, q


def fastq_stats(outdir, paths):
    rows_cycle, rows_trans, rows_runs, rows_perread, rows_n, summary = [], [], [], [], [], []
    all_chars, tops = set(), set()
    for mate, path in enumerate(paths, start=1):
        cyc = []  # cyc[c][char] = count
        trans = []  # trans[c][(a, b)] = count, a at cycle c, b at c+1
        runs = {}  # (char, length) -> count
        lowq = {}  # per-read counts of each quality below the top -> reads
        ncyc = {}  # (cycle, qual char) -> N calls
        nqual = {}  # qual char of N -> count
        alphabet = {}
        reads = bases = n_calls = 0
        maxlen = 0
        for s, q in fastq_records(path):
            reads += 1
            L = len(q)
            bases += L
            if L > maxlen:
                for _ in range(L - maxlen):
                    cyc.append({})
                    trans.append({})
                maxlen = L
            for c, ch in enumerate(q):
                d = cyc[c]
                d[ch] = d.get(ch, 0) + 1
            for c in range(L - 1):
                k = (q[c], q[c + 1])
                d = trans[c]
                d[k] = d.get(k, 0) + 1
            for ch, g in itertools.groupby(q):
                k = (ch, sum(1 for _ in g))
                runs[k] = runs.get(k, 0) + 1
            if "N" in s:
                for c, b in enumerate(s):
                    if b == "N":
                        n_calls += 1
                        k = (c, q[c])
                        ncyc[k] = ncyc.get(k, 0) + 1
                        nqual[q[c]] = nqual.get(q[c], 0) + 1
        # alphabet and top quality
        for d in cyc:
            for ch, n in d.items():
                alphabet[ch] = alphabet.get(ch, 0) + n
        chars = sorted(alphabet)
        top = chars[-1]
        # second pass for per-read low-quality counts (cheap: re-read)
        for _, q in fastq_records(path):
            k = tuple(sorted((ch, q.count(ch)) for ch in set(q) if ch != top))
            lowq[k] = lowq.get(k, 0) + 1
        for c, d in enumerate(cyc):
            tot = sum(d.values())
            rows_cycle.append((mate, c + 1, tot, d))
        for c, d in enumerate(trans):
            rows_trans.append((mate, c + 1, d))
        all_chars.update(chars)
        for (ch, ln), n in sorted(runs.items()):
            rows_runs.append([mate, ch, ord(ch) - 33, ln, n])
        for k, n in lowq.items():
            rows_perread.append((mate, dict(k), n))
        tops.add(top)
        for (c, ch), n in sorted(ncyc.items()):
            rows_n.append([mate, c + 1, ch, ord(ch) - 33, n])
        qsum = sum((ord(ch) - 33) * n for ch, n in alphabet.items())
        q30 = sum(n for ch, n in alphabet.items() if ord(ch) - 33 >= 30)
        summary += [
            [mate, "reads", reads],
            [mate, "bases", bases],
            [mate, "max_length", maxlen],
            [mate, "alphabet", "".join(chars)],
            [mate, "alphabet_phred", ",".join(str(ord(ch) - 33) for ch in chars)],
        ]
        for ch in chars:
            summary.append([mate, f"frac_Q{ord(ch) - 33}", fmt(alphabet[ch] / bases)])
        summary += [
            [mate, "mean_q", fmt(qsum / bases)],
            [mate, "frac_ge_q30", fmt(q30 / bases)],
            [mate, "n_calls", n_calls],
            [mate, "n_rate", fmt(n_calls / bases)],
            [mate, "n_call_quals", ",".join(f"{ch}:{n}" for ch, n in sorted(nqual.items()))],
        ]
    # one column per quality seen in either mate
    chars = sorted(all_chars)
    write_tsv(
        os.path.join(outdir, "cycle_quals.tsv"),
        ["mate", "cycle", "bases"] + [f"Q{ord(ch) - 33}" for ch in chars],
        [[m, c, tot] + [d.get(ch, 0) for ch in chars] for m, c, tot, d in rows_cycle],
    )
    write_tsv(
        os.path.join(outdir, "transitions.tsv"),
        ["mate", "cycle", "from", "to", "count"],
        [[m, c, ord(a) - 33, ord(b) - 33, d.get((a, b), 0)] for m, c, d in rows_trans for a in chars for b in chars],
    )
    write_tsv(os.path.join(outdir, "run_lengths.tsv"), ["mate", "char", "phred", "length", "runs"], rows_runs)
    low = [ch for ch in chars if ch not in tops]
    write_tsv(
        os.path.join(outdir, "read_composition.tsv"),
        ["mate"] + [f"n_Q{ord(ch) - 33}" for ch in low] + ["reads"],
        sorted([m] + [d.get(ch, 0) for ch in low] + [n] for m, d, n in rows_perread),
    )
    write_tsv(os.path.join(outdir, "n_calls.tsv"), ["mate", "cycle", "char", "phred", "count"], rows_n)
    write_tsv(os.path.join(outdir, "fastq_summary.tsv"), ["mate", "key", "value"], summary)


# ---------------------------------------------------------------- sam

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MD_RE = re.compile(r"(\d+)|(\^[A-Z]+)|([A-Z])")


def md_flags(md):
    """Per aligned (M) base: '' for a reference match, else the reference base."""
    out = []
    for num, dele, mis in MD_RE.findall(md):
        if num:
            out.extend([""] * int(num))
        elif mis:
            out.append(mis)
    return out


class Rec:
    __slots__ = ("flag", "rname", "pos", "mapq", "cigar", "seq", "qual", "md", "tlen")

    def __init__(self, f):
        self.flag = int(f[1])
        self.rname = f[2]
        self.pos = int(f[3])
        self.mapq = int(f[4])
        self.cigar = f[5]
        self.tlen = int(f[8])
        self.seq = f[9]
        self.qual = f[10]
        self.md = None
        for t in f[11:]:
            if t.startswith("MD:Z:"):
                self.md = t[5:]

    def mate(self):
        return 1 if self.flag & 0x40 else 2

    def usable(self):
        return not (self.flag & 0x904) and self.mapq >= MIN_MAPQ and self.md is not None


def aligned_bases(rec):
    """Yield (ref_pos, query_index, cycle, base, phred, ref_mismatch_base)."""
    flags = md_flags(rec.md)
    L = len(rec.seq)
    rev = rec.flag & 0x10
    qi, rp, k = 0, rec.pos, 0
    for n, op in CIGAR_RE.findall(rec.cigar):
        n = int(n)
        if op in "M=X":
            for _ in range(n):
                cyc = L - 1 - qi if rev else qi
                yield rp, qi, cyc, rec.seq[qi], ord(rec.qual[qi]) - 33, flags[k]
                qi += 1
                rp += 1
                k += 1
        elif op in "IS":
            qi += n
        elif op in "DN":
            rp += n


def indels(rec):
    """Return (ref_start, ref_end, {(op, ref_pos, length)}) of an alignment."""
    rp, out = rec.pos, set()
    for n, op in CIGAR_RE.findall(rec.cigar):
        n = int(n)
        if op in "M=X":
            rp += n
        elif op == "I":
            out.add(("I", rp, n))
        elif op == "D":
            out.add(("D", rp, n))
            rp += n
        elif op == "N":
            rp += n
    return rec.pos, rp, out


# Indels closer than this to either mate's alignment end are not compared:
# the aligner may soft-clip there instead.
INDEL_MARGIN = 10


def sam_stats(outdir, stream):
    # err[mate][cycle][phred] = [aligned, mismatches]
    err = {1: {}, 2: {}}
    subs = {}  # (phred, reference base, read base), both in read orientation -> count
    nalign = {1: [0, 0], 2: [0, 0]}  # N bases among aligned positions, aligned positions
    indel = {m: {"aligned": 0, "ins_events": 0, "ins_bases": 0, "del_events": 0, "del_bases": 0, "reads": 0} for m in (1, 2)}
    ov_bases, ov_err = {}, {}  # (mate, phred) -> count
    ov = {"positions": 0, "concordant_nonref": 0, "discordant_both_nonref": 0, "pairs": 0,
          "indel_positions": 0, "indel_concordant": 0, "indel_discordant": 0}
    inserts = {}
    counts = {"records": 0, "primary": 0, "primary_mapped": 0, "usable": 0, "pairs": 0, "proper_usable_pairs": 0}

    def account(rec):
        m = rec.mate()
        e = err[m]
        ix = indel[m]
        ix["reads"] += 1
        for n, op in CIGAR_RE.findall(rec.cigar):
            if op == "I":
                ix["ins_events"] += 1
                ix["ins_bases"] += int(n)
            elif op == "D":
                ix["del_events"] += 1
                ix["del_bases"] += int(n)
            elif op in "M=X":
                ix["aligned"] += int(n)
        for _, _, cyc, base, q, refb in aligned_bases(rec):
            nalign[m][1] += 1
            if base == "N":
                nalign[m][0] += 1
                continue
            if refb == "N":
                continue
            d = e.setdefault(cyc, {})
            a = d.setdefault(q, [0, 0])
            a[0] += 1
            ref = refb or base
            if rec.flag & 0x10:
                k = (q, COMP.get(ref, ref), COMP.get(base, base))
            else:
                k = (q, ref, base)
            subs[k] = subs.get(k, 0) + 1
            if refb:
                a[1] += 1

    def overlap(r1, r2):
        pos1 = {}
        for rp, _, _, base, q, refb in aligned_bases(r1):
            pos1[rp] = (base, q, refb)
        hit = False
        for rp, _, _, base2, q2, refb2 in aligned_bases(r2):
            t = pos1.get(rp)
            if t is None:
                continue
            base1, q1, refb1 = t
            if base1 == "N" or base2 == "N" or refb1 == "N" or refb2 == "N":
                continue
            hit = True
            ov["positions"] += 1
            for k in ((1, q1), (2, q2)):
                ov_bases[k] = ov_bases.get(k, 0) + 1
            if refb1 and refb2:
                if base1 == base2:
                    ov["concordant_nonref"] += 1
                else:
                    ov["discordant_both_nonref"] += 1
            elif refb1:
                ov_err[(1, q1)] = ov_err.get((1, q1), 0) + 1
            elif refb2:
                ov_err[(2, q2)] = ov_err.get((2, q2), 0) + 1
        if hit:
            ov["pairs"] += 1
        s1, e1, i1 = indels(r1)
        s2, e2, i2 = indels(r2)
        lo, hi = max(s1, s2) + INDEL_MARGIN, min(e1, e2) - INDEL_MARGIN
        if hi > lo:
            ov["indel_positions"] += hi - lo
            for mine, other in ((i1, i2), (i2, i1)):
                for ev in mine:
                    if lo <= ev[1] < hi:
                        if ev in other:
                            ov["indel_concordant"] += 1  # counted once per mate
                        else:
                            ov["indel_discordant"] += 1

    def pair(group):
        counts["pairs"] += 1
        prim = {}
        for rec in group:
            if not (rec.flag & 0x900):
                prim[rec.mate()] = rec
        r1, r2 = prim.get(1), prim.get(2)
        for r in (r1, r2):
            if r is not None and r.usable():
                account(r)
        if r1 and r2 and r1.usable() and r2.usable() and r1.flag & 0x2 and r1.rname == r2.rname:
            counts["proper_usable_pairs"] += 1
            t = abs(r1.tlen)
            inserts[t] = inserts.get(t, 0) + 1
            overlap(r1, r2)

    group, name = [], None
    for line in stream:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        counts["records"] += 1
        rec = Rec(f)
        if not (rec.flag & 0x900):
            counts["primary"] += 1
            if not (rec.flag & 0x4):
                counts["primary_mapped"] += 1
                if rec.usable():
                    counts["usable"] += 1
        if f[0] != name:
            if group:
                pair(group)
            group, name = [], f[0]
        group.append(rec)
    if group:
        pair(group)

    # per cycle x quality
    rows = []
    by_q = {}
    for m in (1, 2):
        for cyc in sorted(err[m]):
            for q in sorted(err[m][cyc]):
                a, k = err[m][cyc][q]
                rows.append([m, cyc + 1, q, a, k])
                t = by_q.setdefault((m, q), [0, 0])
                t[0] += a
                t[1] += k
    write_tsv(os.path.join(outdir, "errors_cycle.tsv"), ["mate", "cycle", "phred", "aligned", "mismatches"], rows)
    write_tsv(
        os.path.join(outdir, "substitutions.tsv"),
        ["phred", "ref", "aligned", "to_A", "to_C", "to_G", "to_T"],
        [
            [q, ref, sum(subs.get((q, ref, b), 0) for b in "ACGT")] + [subs.get((q, ref, b), 0) if b != ref else 0 for b in "ACGT"]
            for q in sorted({k[0] for k in subs})
            for ref in "ACGT"
        ],
    )

    # Sample variants: the Q37 mismatch rate less the Q37 error rate where
    # mates overlap (variants do not depend on the quality, so this rate is
    # subtracted from every bin). The concordant non-reference rate in the
    # overlap is reported too; it comes from short fragments only.
    pos = ov["positions"]
    concordant = ov["concordant_nonref"] / pos if pos else float("nan")
    a37 = sum(by_q.get((m, 37), [0, 0])[0] for m in (1, 2))
    k37 = sum(by_q.get((m, 37), [0, 0])[1] for m in (1, 2))
    ob37 = sum(ov_bases.get((m, 37), 0) for m in (1, 2))
    oe37 = sum(ov_err.get((m, 37), 0) for m in (1, 2))
    snp = k37 / a37 - oe37 / ob37 if a37 and ob37 else float("nan")
    rows = []
    for (m, q) in sorted(by_q):
        a, k = by_q[(m, q)]
        raw = k / a if a else float("nan")
        ob, oe = ov_bases.get((m, q), 0), ov_err.get((m, q), 0)
        rows.append([m, q, a, k, fmt(raw), fmt(max(raw - snp, 0)), ob, oe, fmt(oe / ob) if ob else "nan", fmt(10 ** (-q / 10))])
    for q in sorted({q for _, q in by_q}):
        a = sum(by_q.get((m, q), [0, 0])[0] for m in (1, 2))
        k = sum(by_q.get((m, q), [0, 0])[1] for m in (1, 2))
        ob = sum(ov_bases.get((m, q), 0) for m in (1, 2))
        oe = sum(ov_err.get((m, q), 0) for m in (1, 2))
        raw = k / a if a else float("nan")
        rows.append(["both", q, a, k, fmt(raw), fmt(max(raw - snp, 0)), ob, oe, fmt(oe / ob) if ob else "nan", fmt(10 ** (-q / 10))])
    write_tsv(
        os.path.join(outdir, "errors_quality.tsv"),
        ["mate", "phred", "aligned", "mismatches", "raw_rate", "rate_minus_variants", "overlap_bases", "overlap_errors", "overlap_rate", "phred_rate"],
        rows,
    )

    # insert sizes
    write_tsv(os.path.join(outdir, "insert_sizes.tsv"), ["insert", "pairs"], sorted(inserts.items()))
    vals = []
    for t, n in inserts.items():
        vals.extend([t] * n)
    vals.sort()

    def pct(p):
        return vals[min(len(vals) - 1, int(p * len(vals)))] if vals else 0

    summary = [[k, v] for k, v in counts.items()]
    if vals:
        mu = statistics.fmean(vals)
        sd = statistics.pstdev(vals)
        # robust spread: MAD-based SD and the central 98 % (as in Picard's
        # CollectInsertSizeMetrics, which trims to median +/- 10 MAD)
        med = statistics.median(vals)
        mad = statistics.median([abs(v - med) for v in vals])
        core = [v for v in vals if abs(v - med) <= 10 * mad]
        lv = [math.log(v) for v in core if v > 0]
        lm, ls = statistics.fmean(lv), statistics.pstdev(lv)
        cm = statistics.fmean(core)
        skew = sum((v - cm) ** 3 for v in core) / len(core) / statistics.pstdev(core) ** 3
        summary += [
            ["insert_pairs", len(vals)],
            ["insert_mean", fmt(mu)],
            ["insert_sd", fmt(sd)],
            ["insert_median", med],
            ["insert_mad", mad],
            ["insert_p01", pct(0.01)],
            ["insert_p05", pct(0.05)],
            ["insert_p25", pct(0.25)],
            ["insert_p75", pct(0.75)],
            ["insert_p95", pct(0.95)],
            ["insert_p99", pct(0.99)],
            ["insert_core_pairs", len(core)],
            ["insert_core_mean", fmt(cm)],
            ["insert_core_sd", fmt(statistics.pstdev(core))],
            ["insert_core_skewness", fmt(skew)],
            ["insert_core_log_mean", fmt(lm)],
            ["insert_core_log_sd", fmt(ls)],
            ["insert_below_151", sum(1 for v in vals if v < 151)],
            ["insert_below_302", sum(1 for v in vals if v < 302)],
        ]
    for m in (1, 2):
        ix = indel[m]
        a = ix["aligned"]
        summary += [
            [f"mate{m}_aligned_bases", a],
            [f"mate{m}_ins_events_per_base", fmt(ix["ins_events"] / a)],
            [f"mate{m}_del_events_per_base", fmt(ix["del_events"] / a)],
            [f"mate{m}_ins_bases_per_base", fmt(ix["ins_bases"] / a)],
            [f"mate{m}_del_bases_per_base", fmt(ix["del_bases"] / a)],
            [f"mate{m}_aligned_n_rate", fmt(nalign[m][0] / nalign[m][1])],
        ]
        tot = sum(v[0] for (mm, _), v in by_q.items() if mm == m)
        mis = sum(v[1] for (mm, _), v in by_q.items() if mm == m)
        summary.append([f"mate{m}_raw_mismatch_rate", fmt(mis / tot)])
    summary += [
        ["overlap_pairs", ov["pairs"]],
        ["overlap_positions", pos],
        ["overlap_concordant_nonref", ov["concordant_nonref"]],
        ["overlap_discordant_both_nonref", ov["discordant_both_nonref"]],
        ["concordant_nonref_rate_in_overlap", fmt(concordant)],
        ["variant_rate", fmt(snp)],
        ["overlap_indel_positions", ov["indel_positions"]],
        ["overlap_indel_concordant_events", ov["indel_concordant"] // 2],
        ["overlap_indel_discordant_events", ov["indel_discordant"]],
        ["indel_error_rate_from_overlap", fmt(ov["indel_discordant"] / (2 * ov["indel_positions"])) if ov["indel_positions"] else "nan"],
        ["indel_variant_rate_from_overlap", fmt(ov["indel_concordant"] / 2 / ov["indel_positions"]) if ov["indel_positions"] else "nan"],
    ]
    write_tsv(os.path.join(outdir, "sam_summary.tsv"), ["key", "value"], summary)


def main(argv):
    if len(argv) < 3:
        raise SystemExit(__doc__)
    cmd, outdir = argv[1], argv[2]
    os.makedirs(outdir, exist_ok=True)
    if cmd == "fastq" and len(argv) == 5:
        fastq_stats(outdir, argv[3:5])
    elif cmd == "sam" and len(argv) == 3:
        sam_stats(outdir, sys.stdin)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv)
