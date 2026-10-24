"""Measure Illumina quality, error and insert-size statistics (#27).

Run with `python3 -I` (no third-party packages). Subcommands:

  fastq OUTDIR R1.fastq R2.fastq
      Quality alphabet, per-cycle quality frequencies, per-cycle transition
      counts P(q[c+1] | q[c]), run lengths, low-quality bases per read and
      N calls, for each mate.

  context OUTDIR R1.fastq R2.fastq FIT.tsv
      Quality frequencies by called-base sequence context (the two bases
      before, the base and the base after) and 10-cycle block, and the
      Q11 and Q25 calls by context within each read quality class of
      FIT.tsv (written by fit_illumina.py), against their expected numbers
      without a context effect (#35).

  pairs OUTDIR R1.fastq R2.fastq FIT.tsv
      How the two mates' read classes go together: the joint distribution
      of read 1's and read 2's quality class (FIT.tsv's per-mate classes),
      fitted by EM to the pairs' Q11 and Q25 counts with the per-mate class
      levels fixed (pair_classes.tsv), and the measured coupling of the
      mates' qualities (pair_summary.tsv, #35).

  structure OUTDIR R1.fastq R2.fastq
      How the qualities are arranged along the reads (#36): transitions of
      order 1 and 2, by cycle block and by read stratum (the reads' count
      of Q11 and Q25 bases), run and segment lengths, isolated dips, the
      kinds of quality change, the excess of low calls at lag k within a
      stratum, each stratum's per-cycle profile and the spread of the
      classifier's block-mean features (structure.tsv).

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


# Cycle blocks of the sequence-context table: 10 cycles each (1-10, 11-20,
# ..., 141-150, 151).
CONTEXT_BLOCK = 10


def read_classes(fit_path):
    """Per mate, the read quality classes of fit.tsv: weights, Q11 and Q25
    multipliers and the expected Q11 and Q25 bases per read."""
    with open(fit_path, encoding="ascii") as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f]
    out = {}
    for mate in ("1", "2"):
        get = lambda p: {int(r["index"]): float(r["value"]) for r in rows if r["param"] == p and r["mate"] == mate}
        w, m11, m25 = get("class_weight"), get("class_mult_q11"), get("class_mult_q25")
        lam = (get("class_lambda_q11")[0], get("class_lambda_q25")[0])
        out[int(mate)] = [(w[k], m11[k], m25[k]) for k in sorted(w)], lam
    return out


def class_posterior(classes, lam, n11, n25):
    """The posterior probabilities of the read classes for a read with n11
    Q11 and n25 Q25 bases under the Poisson mixture of fit_illumina.py, and
    the most probable class."""
    lp = []
    for w, m11, m25 in classes:
        a, b = m11 * lam[0], m25 * lam[1]
        lp.append(math.log(w) + n11 * math.log(a) - a + n25 * math.log(b) - b)
    mx = max(lp)
    e = [math.exp(v - mx) for v in lp]
    z = sum(e)
    return [v / z for v in e], lp.index(mx)


# Read GC bins of the class-by-read table: floor(20 * GC), GC the fraction
# of G and C among the read's A, C, G and T bases; longest homopolymer
# (run of one of A, C, G, T) capped at HP_CAP.
GC_BINS = 20
HP_CAP = 20
HP_RE = re.compile(r"A+|C+|G+|T+")


def longest_homopolymer(s):
    return max((m.end() - m.start() for m in HP_RE.finditer(s)), default=0)


def context_stats(outdir, paths, fit_path):
    """Quality bins by called-base context (#35).

    context_quals.tsv: one row per mate, cycle block and context (the called
    bases two before, one before, at and one after the cycle; '-' past the
    read's ends), with the count of each quality value. Every marginal
    (called base, preceding base, GG motifs, by block) is a sum of rows.

    context_classes.tsv: one row per mate, read class and context: the
    bases, the Q11 and Q25 calls, and their expected numbers if the read's
    calls fell on its cycles in proportion to the mate's per-cycle Q11 and
    Q25 fractions whatever the context (exp = sum over the read's cycles in
    that context of n_b f_b(c) / sum_c f_b(c), n_b the read's calls in bin
    b). A read's class is its most probable class (class_posterior) given
    its Q11 and Q25 counts. obs/exp is the context's factor within reads of a
    class; conditioning on the read's own counts removes the class's level.

    class_reads.tsv: per mate, read GC bin (floor(20 GC)) and longest
    homopolymer (capped at HP_CAP), the reads and the sum of their class
    posterior probabilities: P(class | GC, homopolymer).
    """
    fits = read_classes(fit_path)
    marg = {}  # (mate, block, ctx4, qual char) -> count
    cls = {}  # (mate, class, ctx4) -> [bases, n11, n25, exp11, exp25]
    gcs = {}  # (mate, GC bin, homopolymer) -> [reads, posterior sum of each class]
    chars = set()
    maxlen = 0
    for mate, path in enumerate(paths, start=1):
        classes, lam = fits[mate]
        frac = []  # per cycle: [Q11, Q25] calls
        for s, q in fastq_records(path):
            while len(frac) < len(q):
                frac.append([0, 0])
            for c, ch in enumerate(q):
                if ch == ",":
                    frac[c][0] += 1
                elif ch == ":":
                    frac[c][1] += 1
        cache = {}
        for s, q in fastq_records(path):
            L = len(q)
            maxlen = max(maxlen, L)
            n11, n25 = q.count(","), q.count(":")
            key = (n11, n25)
            if key not in cache:
                cache[key] = class_posterior(classes, lam, n11, n25)
            post, k = cache[key]
            acgt = sum(s.count(b) for b in "ACGT")
            if acgt:
                gk = (mate, (s.count("G") + s.count("C")) * GC_BINS // acgt, min(longest_homopolymer(s), HP_CAP))
                g = gcs.get(gk)
                if g is None:
                    g = gcs[gk] = [0] + [0.0] * len(classes)
                g[0] += 1
                for i, v in enumerate(post):
                    g[1 + i] += v
            t11 = sum(f[0] for f in frac[:L]) or 1
            t25 = sum(f[1] for f in frac[:L]) or 1
            u11, u25 = n11 / t11, n25 / t25
            p = "--" + s + "-"
            for c, ch in enumerate(q):
                ctx = p[c : c + 4]
                mk = (mate, c // CONTEXT_BLOCK, ctx, ch)
                marg[mk] = marg.get(mk, 0) + 1
                ck = (mate, k, ctx)
                v = cls.get(ck)
                if v is None:
                    v = cls[ck] = [0, 0, 0, 0.0, 0.0]
                v[0] += 1
                if ch == ",":
                    v[1] += 1
                elif ch == ":":
                    v[2] += 1
                v[3] += u11 * frac[c][0]
                v[4] += u25 * frac[c][1]
        chars.update(mk[3] for mk in marg)
    chars = sorted(chars)
    cells = {}
    for (mate, block, ctx, ch), n in marg.items():
        cells.setdefault((mate, block, ctx), {})[ch] = n
    write_tsv(
        os.path.join(outdir, "context_quals.tsv"),
        ["mate", "first_cycle", "last_cycle", "prev2", "prev", "base", "next"] + [f"Q{ord(ch) - 33}" for ch in chars],
        [
            [mate, block * CONTEXT_BLOCK + 1, min((block + 1) * CONTEXT_BLOCK, maxlen), ctx[0], ctx[1], ctx[2], ctx[3]]
            + [d.get(ch, 0) for ch in chars]
            for (mate, block, ctx), d in sorted(cells.items())
        ],
    )
    write_tsv(
        os.path.join(outdir, "context_classes.tsv"),
        ["mate", "class", "prev2", "prev", "base", "next", "bases", "n_Q11", "n_Q25", "exp_Q11", "exp_Q25"],
        [[m, k, x[0], x[1], x[2], x[3], v[0], v[1], v[2], fmt(v[3]), fmt(v[4])] for (m, k, x), v in sorted(cls.items())],
    )
    nk = max(len(v) for v in gcs.values()) - 1
    write_tsv(
        os.path.join(outdir, "class_reads.tsv"),
        ["mate", "gc_bin", "gc_from", "homopolymer", "reads"] + [f"class{k}" for k in range(nk)],
        [[m, b, fmt(b / GC_BINS), h, v[0]] + [fmt(x) for x in v[1:]] for (m, b, h), v in sorted(gcs.items())],
    )


# ---------------------------------------------------------------- pairs

PAIR_EM_ITERS = 300


def pair_stats(outdir, paths, fit_path):
    """Joint read classes of the mates (#35).

    The classes of each mate (fit_illumina.py, fitted to that mate's reads
    alone) give every read a likelihood under each class: Poisson Q11 and
    Q25 counts with means m_kb * lambda_b. The pair's classes (k1, k2) get
    weights W[k1][k2], fitted by EM to the pairs' four counts with those
    likelihoods fixed, starting from independent mates (w1[k1] w2[k2]).
    pair_classes.tsv holds W (and the log-likelihoods of independent and
    coupled mates); fit_illumina.py rakes it to the per-mate weights.
    pair_summary.tsv holds the measured coupling: the correlation of the
    mates' Q11 + Q25 + Q2 counts and of their Q11 + Q2 counts (Phred < 20).
    """
    fits = read_classes(fit_path)
    pairs = {}
    low30, low20 = ([], []), ([], [])
    for (s1, q1), (s2, q2) in zip(fastq_records(paths[0]), fastq_records(paths[1])):
        k = (q1.count(","), q1.count(":"), q2.count(","), q2.count(":"))
        pairs[k] = pairs.get(k, 0) + 1
        for i, q in enumerate((q1, q2)):
            n2 = q.count("#")
            low30[i].append(k[2 * i] + k[2 * i + 1] + n2)
            low20[i].append(k[2 * i] + n2)
    rows = []
    for (a, b, c, d), n in sorted(pairs.items()):
        ls = []
        for mate, (x, y) in ((1, (a, b)), (2, (c, d))):
            classes, lam = fits[mate]
            l = []
            for _, m11, m25 in classes:
                u, v = m11 * lam[0], m25 * lam[1]
                l.append(x * math.log(u) - u + y * math.log(v) - v - math.lgamma(x + 1) - math.lgamma(y + 1))
            mx = max(l)
            ls.append((mx, [math.exp(t - mx) for t in l]))
        rows.append((ls, n))
    w1 = [w for w, _, _ in fits[1][0]]
    w2 = [w for w, _, _ in fits[2][0]]
    K1, K2 = len(w1), len(w2)
    W = [[w1[i] * w2[j] for j in range(K2)] for i in range(K1)]
    tot = sum(n for _, n in rows)

    def step(W):
        S = [[0.0] * K2 for _ in range(K1)]
        ll = 0.0
        for ((m1, e1), (m2, e2)), n in rows:
            p = [[W[i][j] * e1[i] * e2[j] for j in range(K2)] for i in range(K1)]
            z = sum(map(sum, p))
            ll += n * (m1 + m2 + math.log(z))
            for i in range(K1):
                for j in range(K2):
                    S[i][j] += n * p[i][j] / z
        return [[S[i][j] / tot for j in range(K2)] for i in range(K1)], ll

    _, ll_indep = step(W)
    for _ in range(PAIR_EM_ITERS):
        W, ll = step(W)
    write_tsv(
        os.path.join(outdir, "pair_classes.tsv"),
        ["class_r1", "class_r2", "weight"],
        [[i, j, fmt(W[i][j], 8)] for i in range(K1) for j in range(K2)],
    )
    write_tsv(
        os.path.join(outdir, "pair_summary.tsv"),
        ["key", "value"],
        [
            ["pairs", tot],
            ["distinct_count_tuples", len(rows)],
            ["loglik_independent", fmt(ll_indep, 10)],
            ["loglik_coupled", fmt(ll, 10)],
            ["em_iterations", PAIR_EM_ITERS],
            ["lowq30_count_correlation", fmt(statistics.correlation(*low30))],
            ["lowq20_count_correlation", fmt(statistics.correlation(*low20))],
        ],
    )


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


# ---------------------------------------------------------------- structure

# Reads are grouped by their number of Q11 and Q25 bases (STRATA, [lo, hi)),
# a stand-in for the read class that both real and generated reads have.
STRATA = [(0, 1), (1, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 10**6)]
LAGS = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]
# Run (and non-Q37 segment) length bins, [lo, hi).
RUN_BINS = [1, 2, 3, 4, 5, 10, 20, 40, 80, 10**6]
# Cycle blocks of the per-stratum low-call profile (0-based [lo, hi); the
# last two are the cycles before the final one and the final cycle).
SHAPE_BLOCKS = [(0, 1), (1, 3), (3, 10), (10, 30), (30, 60), (60, 90), (90, 120), (120, 140)]
TRANS_BLOCK = 30  # cycle blocks of the transition table
QUAL_BLOCKS = 10  # the classifier's equal cycle blocks (internal/testing)
BINS = "#,:F"  # Q2, Q11, Q25, Q37


def stratum(n):
    for i, (lo, hi) in enumerate(STRATA):
        if lo <= n < hi:
            return i
    return len(STRATA) - 1


def bin_label(lo, hi):
    return f"{lo}-{hi - 1}" if hi < 10**6 else f">={lo}"


def structure_stats(outdir, paths):
    """How the qualities of a read set are arranged along the read (#36).

    Writes structure.tsv, one row per statistic: section, mate, key, value.
    Sections:
      trans1      P(q[c+1] | q[c]) pooled over cycles (key from>to);
      trans2      P(q[c+2] | q[c], q[c+1]) (key q0q1>q2), for q0, q1 below Q37
                  or Q37 (second order: is the chain order 1?);
      trans_block P(q[c+1] | q[c]) per block of TRANS_BLOCK cycles (key
                  from>to@first-last), low bins only;
      trans_stratum P(q[c+1] | q[c]) within the reads of each stratum (reads
                  with n Q11 + Q25 bases, STRATA; key from>to@n), low bins;
      runs        maximal runs of one bin per read, by length (key bin@len);
      segments    maximal runs of bins below Q37 per read, by length, the
                  share of them ending at the last cycle and their median
                  start cycle;
      dips        isolated low calls per read: Q37, x, Q37 (key x);
      changes     quality changes per read by type (key from>to), and the
                  mean fraction of adjacent cycles that differ
                  (qual_change_frac, the classifier's feature);
      changes_per_low  changes per Q11 or Q25 base within each stratum;
      lag_excess  within each stratum, P(low at c+k | low at c) / P(low at
                  c+k) for low = Q11 or Q25 (key n@k): 1 for cycles that are
                  independent given the stratum;
      shape       within each stratum, the low-call fraction of cycle blocks
                  relative to the stratum's mean (key n@first-last): a
                  common shape for every stratum means the read classes
                  share one per-cycle curve;
      strata      the share of reads in each stratum;
      block_sd    the standard deviation over reads of the mean quality of
                  each of QUAL_BLOCKS equal cycle blocks, and of the
                  longest and the final run below Q20 (classifier features).
    """
    rows = []
    for mate, path in enumerate(paths, start=1):
        reads = 0
        t1 = {}
        t2 = {}
        tb = {}
        ts = {}
        runs = {}
        seg_len = {}
        seg_end = {}
        seg_start = {}
        dips = {",": 0, ":": 0}
        changes = {}
        change_frac = 0.0
        ch_low = {}  # stratum -> [changes, lows]
        nstrat = [0] * len(STRATA)
        # lag excess: per stratum, per lag: [pairs both low, lows at c with c+k < L]
        lag_pairs = {}
        lag_from = {}
        low_by_cycle = {}  # stratum -> per-cycle low counts
        bsum = [0.0] * QUAL_BLOCKS
        bsq = [0.0] * QUAL_BLOCKS
        longest = [0.0, 0.0]
        tail = [0.0, 0.0]
        maxlen = 0
        for _, q in fastq_records(path):
            reads += 1
            L = len(q)
            maxlen = max(maxlen, L)
            lows = [c for c, ch in enumerate(q) if ch == "," or ch == ":"]
            n = len(lows)
            s = stratum(n)
            nstrat[s] += 1
            lc = low_by_cycle.setdefault(s, [])
            while len(lc) < L:
                lc.append(0)
            for c in lows:
                lc[c] += 1
            for a, b in zip(q, q[1:]):
                k = a + b
                t1[k] = t1.get(k, 0) + 1
            for c in range(L - 1):
                if q[c] in ",:":
                    k = (c // TRANS_BLOCK, q[c], q[c + 1])
                    tb[k] = tb.get(k, 0) + 1
                    k = (s, q[c], q[c + 1])
                    ts[k] = ts.get(k, 0) + 1
            for a, b, c3 in zip(q, q[1:], q[2:]):
                k = a + b + c3
                t2[k] = t2.get(k, 0) + 1
            nch = 0
            for ch, g in itertools.groupby(q):
                k = (ch, sum(1 for _ in g))
                runs[k] = runs.get(k, 0) + 1
            c = 0
            while c < L:
                if q[c] != "F":
                    e = c
                    while e < L and q[e] != "F":
                        e += 1
                    ln = e - c
                    seg_len[ln] = seg_len.get(ln, 0) + 1
                    if e == L:
                        seg_end[ln] = seg_end.get(ln, 0) + 1
                    seg_start.setdefault(ln, []).append(c)
                    c = e
                else:
                    c += 1
            for a, b, c3 in zip(q, q[1:], q[2:]):
                if a == "F" and c3 == "F" and b in dips:
                    dips[b] += 1
            for a, b in zip(q, q[1:]):
                if a != b:
                    nch += 1
                    k = a + b
                    changes[k] = changes.get(k, 0) + 1
            if L > 1:
                change_frac += nch / (L - 1)
            v = ch_low.setdefault(s, [0, 0])
            v[0] += nch
            v[1] += n
            ls = set(lows)
            for k in LAGS:
                kk = (s, k)
                lag_from[kk] = lag_from.get(kk, 0) + sum(1 for c in lows if c + k < L)
                lag_pairs[kk] = lag_pairs.get(kk, 0) + sum(1 for c in lows if c + k in ls)
            for b in range(QUAL_BLOCKS):
                lo, hi = b * L // QUAL_BLOCKS, (b + 1) * L // QUAL_BLOCKS
                m = sum(ord(ch) - 33 for ch in q[lo:hi]) / max(hi - lo, 1)
                bsum[b] += m
                bsq[b] += m * m
            run = best = 0
            for ch in q:
                if ord(ch) - 33 < 20:
                    run += 1
                    best = max(best, run)
                else:
                    run = 0
            longest[0] += best
            longest[1] += best * best
            tail[0] += run
            tail[1] += run * run
        m = str(mate)

        def add(section, key, value):
            rows.append([section, m, key, fmt(value)])

        for a in BINS:
            tot = sum(t1.get(a + b, 0) for b in BINS)
            for b in BINS:
                if tot:
                    add("trans1", f"Q{ord(a) - 33}>Q{ord(b) - 33}", t1.get(a + b, 0) / tot)
        for a in ",:F":
            for b in ",:F":
                tot = sum(t2.get(a + b + x, 0) for x in BINS)
                for x in ",:F":
                    if tot:
                        add("trans2", f"Q{ord(a) - 33}Q{ord(b) - 33}>Q{ord(x) - 33}", t2.get(a + b + x, 0) / tot)
        for blk in sorted({k[0] for k in tb}):
            lo, hi = blk * TRANS_BLOCK + 1, min((blk + 1) * TRANS_BLOCK, maxlen - 1)
            for a in ",:":
                tot = sum(tb.get((blk, a, b), 0) for b in BINS)
                for b in ",:F":
                    if tot:
                        add("trans_block", f"Q{ord(a) - 33}>Q{ord(b) - 33}@{lo}-{hi}", tb.get((blk, a, b), 0) / tot)
        for si, (lo, hi) in enumerate(STRATA):
            for a in ",:":
                tot = sum(ts.get((si, a, b), 0) for b in BINS)
                for b in ",:F":
                    if tot:
                        add("trans_stratum", f"Q{ord(a) - 33}>Q{ord(b) - 33}@{bin_label(lo, hi)}", ts.get((si, a, b), 0) / tot)
        for ch in BINS:
            for lo, hi in zip(RUN_BINS, RUN_BINS[1:]):
                nr = sum(v for (c2, ln), v in runs.items() if c2 == ch and lo <= ln < hi)
                add("runs", f"Q{ord(ch) - 33}@{bin_label(lo, hi)}", nr / reads)
        for lo, hi in zip(RUN_BINS, RUN_BINS[1:]):
            ns = sum(v for ln, v in seg_len.items() if lo <= ln < hi)
            add("segments", f"per_read@{bin_label(lo, hi)}", ns / reads)
            if ns:
                add("segments", f"ending_at_last_cycle@{bin_label(lo, hi)}", sum(v for ln, v in seg_end.items() if lo <= ln < hi) / ns)
                starts = sorted(x for ln, xs in seg_start.items() if lo <= ln < hi for x in xs)
                add("segments", f"median_start_cycle@{bin_label(lo, hi)}", starts[len(starts) // 2] + 1)
        add("segments", "per_read", sum(seg_len.values()) / reads)
        for ch, v in dips.items():
            add("dips", f"Q{ord(ch) - 33}", v / reads)
        for a in BINS:
            for b in BINS:
                if a != b:
                    add("changes", f"Q{ord(a) - 33}>Q{ord(b) - 33}", changes.get(a + b, 0) / reads)
        add("changes", "qual_change_frac", change_frac / reads)
        for si, (lo, hi) in enumerate(STRATA):
            v = ch_low.get(si)
            if v and v[1]:
                add("changes_per_low", bin_label(lo, hi), v[0] / v[1])
            add("strata", bin_label(lo, hi), nstrat[si] / reads)
            lc = low_by_cycle.get(si)
            if not lc or not sum(lc):
                continue
            L = len(lc)
            mean = sum(lc) / L
            for k in LAGS:
                # P(low at c+k | low at c) / P(low at c', c' >= k)
                frm = lag_from.get((si, k), 0)
                later = sum(lc[k:]) / (nstrat[si] * (L - k))
                if frm and later:
                    add("lag_excess", f"{bin_label(lo, hi)}@{k}", lag_pairs[(si, k)] / frm / later)
            blocks = SHAPE_BLOCKS + [(SHAPE_BLOCKS[-1][1], L - 1), (L - 1, L)]
            for a, b in blocks:
                if b > a:
                    add("shape", f"{bin_label(lo, hi)}@{a + 1}-{b}", sum(lc[a:b]) / (b - a) / mean)
        for b in range(QUAL_BLOCKS):
            mu = bsum[b] / reads
            add("block_sd", f"mean_qual_block_{b + 1}/{QUAL_BLOCKS}", math.sqrt(max(bsq[b] / reads - mu * mu, 0)))
        for name, v in (("lowq_longest_run", longest), ("lowq_tail_run", tail)):
            mu = v[0] / reads
            add("block_sd", f"{name}_mean", mu)
            add("block_sd", f"{name}_sd", math.sqrt(max(v[1] / reads - mu * mu, 0)))
    write_tsv(os.path.join(outdir, "structure.tsv"), ["section", "mate", "key", "value"], rows)


def main(argv):
    if len(argv) < 3:
        raise SystemExit(__doc__)
    cmd, outdir = argv[1], argv[2]
    os.makedirs(outdir, exist_ok=True)
    if cmd == "fastq" and len(argv) == 5:
        fastq_stats(outdir, argv[3:5])
    elif cmd == "context" and len(argv) == 6:
        context_stats(outdir, argv[3:5], argv[5])
    elif cmd == "pairs" and len(argv) == 6:
        pair_stats(outdir, argv[3:5], argv[5])
    elif cmd == "sam" and len(argv) == 3:
        sam_stats(outdir, sys.stdin)
    elif cmd == "structure" and len(argv) == 5:
        structure_stats(outdir, argv[3:5])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv)
