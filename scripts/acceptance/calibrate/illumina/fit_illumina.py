"""Fit the phantom-science Illumina model to the tables of illumina_stats.py (#27).

Run with `python3 -I fit_illumina.py RESULTS_DIR [GO_FILE]`; writes
RESULTS_DIR/fit.tsv (parameter, mate, index, value) and, with GO_FILE, the
parameters as Go source: the run's Profile
(internal/genomics/illumina/profile_<accession>.go), described by
RESULTS_DIR/profile.tsv (name, accession, study, instrument, center,
library, pairs, results; for a profile pooled from several runs also id,
the Go name's suffix, pairs_per_run, and the runs in
RESULTS_DIR/profile_runs.tsv) and the headline statistics of the measured
reads.

Model (see the package documentation of internal/genomics/illumina):

* Per-cycle fraction of each low bin (Q11, Q25), l_b(c), is piecewise linear
  in the cycle between fixed knots (least squares over cycles 1..L-1 of the
  measured reads), extended linearly past the last knot; the final cycle of
  a read gets its own factor (the measured final-cycle fraction over the
  curve's value there).
* Reads differ in quality: each read draws a class k (weight w_k) that
  multiplies the low-bin fractions by m_kb. The classes are fitted by EM to
  the per-read counts of Q11 and Q25 bases, treating them as independent
  Poisson counts given the class (mean m_kb * sum_c l_b(c)); EM keeps
  sum_k w_k m_kb = 1, so the class mixture has the per-cycle fractions l_b(c).
  Given its class, a read's bins are drawn independently per cycle (the
  measured autocorrelation of low-quality calls is flat in the lag, so it
  comes from the read, not from neighbouring cycles).
* Sequence context (#35): within the reads of class k, P(Q11) and P(Q25) at
  a cycle are multiplied by F_kb(x), x the 4-base context (called bases two
  before, one before, at and one after the cycle). F is the context's calls
  over their number expected from the read's own counts and the per-cycle
  fractions (context_classes.tsv), shrunk towards the 3-base and 1-base
  contexts (CONTEXT_PRIOR) and normalised to average 1 over the sample's
  context frequencies. Without RESULTS_DIR/context_classes.tsv only the
  classes are fitted (the first pass of illumina_calibrate.sh).
* Mate coupling (#35): the two reads of a pair draw their classes together,
  (k1, k2) with weight W[k1][k2]. W is the EM fit of illumina_stats.py
  pairs (pair_classes.tsv, the pairs' four counts with the per-mate class
  levels fixed), raked (iterative proportional fitting) so that its
  margins are the per-mate weights exactly: each mate's reads keep their
  class mixture. Without pair_classes.tsv the mates are independent.
* Substitution rate per bin and cycle: the mismatch rate against hg38 minus
  the rate of sample variants. The variant rate is the Q37 mismatch rate
  minus the Q37 error rate measured where mates overlap (a mismatch in one
  mate only). Q37 has one rate (its mismatch rate does not change with
  cycle); Q11 and Q25 follow e(c) = A + B exp(-(c-1)/tau), fitted by
  weighted least squares with tau on a grid (TAUS; a negative tau is a
  convex rise, A then the asymptote towards early cycles).
"""

import math
import os
import sys

KNOTS = [1, 2, 3, 4, 6, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 135, 150]
LOW = (11, 25)
CLASSES = 6
MAX_LOW = 0.98  # cap on P(Q11) + P(Q25) of a class at a cycle (illumina.maxLow)


def read_tsv(path):
    with open(path, encoding="ascii") as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f]


def solve(a, b):
    """Solve the square system a x = b (Gaussian elimination, partial pivoting)."""
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(m[r][i]))
        m[i], m[p] = m[p], m[i]
        for r in range(i + 1, n):
            f = m[r][i] / m[i][i]
            for c in range(i, n + 1):
                m[r][c] -= f * m[i][c]
    x = [0.0] * n
    for i in reversed(range(n)):
        x[i] = (m[i][n] - sum(m[i][c] * x[c] for c in range(i + 1, n))) / m[i][i]
    return x


def hat(c, knots):
    """Piecewise-linear basis values at cycle c (linear extension past the ends)."""
    n = len(knots)
    v = [0.0] * n
    if c <= knots[0]:
        v[0] = 1.0
        return v
    for j in range(n - 1):
        if c <= knots[j + 1] or j == n - 2:
            t = (c - knots[j]) / (knots[j + 1] - knots[j])
            v[j], v[j + 1] = 1 - t, t
            return v
    return v


def interp(c, knots, vals):
    return sum(h * y for h, y in zip(hat(c, knots), vals))


def fit_curve(cycles, fracs, weights):
    n = len(KNOTS)
    ata = [[0.0] * n for _ in range(n)]
    atb = [0.0] * n
    for c, y, w in zip(cycles, fracs, weights):
        h = hat(c, KNOTS)
        for i in range(n):
            if h[i]:
                atb[i] += w * h[i] * y
                for j in range(n):
                    ata[i][j] += w * h[i] * h[j]
    return solve(ata, atb)


def pois_logpmf(k, lam):
    return k * math.log(lam) - lam - math.lgamma(k + 1)


def em(rows, lam, k_classes, iters=400):
    """Poisson-mixture EM over rows of (counts tuple, reads)."""
    tot = sum(n for _, n in rows)
    # initialise multipliers on a log grid of overall quality
    w = [1.0 / k_classes] * k_classes
    m = [[0.2 * 100 ** (i / (k_classes - 1)) for _ in lam] for i in range(k_classes)]
    ll = None
    for _ in range(iters):
        sw = [0.0] * k_classes
        sn = [[0.0] * len(lam) for _ in range(k_classes)]
        ll = 0.0
        for counts, n in rows:
            lp = [math.log(w[k]) + sum(pois_logpmf(x, m[k][b] * lam[b]) for b, x in enumerate(counts)) for k in range(k_classes)]
            mx = max(lp)
            s = sum(math.exp(v - mx) for v in lp)
            ll += n * (mx + math.log(s))
            for k in range(k_classes):
                r = n * math.exp(lp[k] - mx) / s
                sw[k] += r
                for b, x in enumerate(counts):
                    sn[k][b] += r * x
        w = [s / tot for s in sw]
        m = [[sn[k][b] / (sw[k] * lam[b]) if sw[k] else m[k][b] for b in range(len(lam))] for k in range(k_classes)]
    order = sorted(range(k_classes), key=lambda k: sum(m[k]))
    return [w[k] for k in order], [m[k] for k in order], ll


def class_target(m_k, curves, c):
    """Per-class low-bin probabilities at cycle index c, capped as in the Go model."""
    a, b = m_k[0] * curves[0][c], m_k[1] * curves[1][c]
    if a + b > MAX_LOW:
        a, b = a * MAX_LOW / (a + b), b * MAX_LOW / (a + b)
    return a, b


def self_transition(w, m, curves, rho):
    """Model P(b at c+1 | b at c) for the low bins, pooled over cycles.

    With stickiness rho_b a cycle keeps a low bin b with probability rho_b
    and otherwise draws from nu(b) = pi(b) (1 - rho_b) / (1 - sum rho pi),
    which keeps the per-cycle distribution pi.
    """
    L = len(curves[0])
    num, den = [0.0, 0.0], [0.0, 0.0]
    for k, wk in enumerate(w):
        for c in range(L - 1):
            p0 = class_target(m[k], curves, c)
            p1 = class_target(m[k], curves, c + 1)
            z = 1 - rho[0] * p1[0] - rho[1] * p1[1]
            for i in range(2):
                nu = p1[i] * (1 - rho[i]) / z
                num[i] += wk * p0[i] * (rho[i] + (1 - rho[i]) * nu)
                den[i] += wk * p0[i]
    return [num[i] / den[i] for i in range(2)]


def fit_stickiness(w, m, curves, real):
    rho = [0.0, 0.0]
    for _ in range(4):
        for i in range(2):
            lo, hi = 0.0, 0.9
            for _ in range(40):
                mid = (lo + hi) / 2
                r = rho[:]
                r[i] = mid
                if self_transition(w, m, curves, r)[i] < real[i]:
                    lo = mid
                else:
                    hi = mid
            rho[i] = (lo + hi) / 2
    return rho


# The grid of tau of fit_error; a negative tau is a rate that grows ever
# faster along the read (convex), as Q25 does in some runs (#35).
TAUS = [2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300] + [-t for t in (20, 30, 40, 50, 60, 80, 100, 150, 200, 300)]


def fit_error(rows_by_cycle):
    """Fit e(c) = A + B exp(-(c-1)/tau) to (cycle, aligned, rate) rows."""
    best = None
    for tau in TAUS:
        ata = [[0.0, 0.0], [0.0, 0.0]]
        atb = [0.0, 0.0]
        for c, n, y in rows_by_cycle:
            x = [1.0, math.exp(-(c - 1) / tau)]
            for i in range(2):
                atb[i] += n * x[i] * y
                for j in range(2):
                    ata[i][j] += n * x[i] * x[j]
        a, b = solve(ata, atb)
        # binomial deviance of the fit
        dev = 0.0
        for c, n, y in rows_by_cycle:
            p = min(max(a + b * math.exp(-(c - 1) / tau), 1e-9), 1 - 1e-9)
            k = y * n
            dev -= k * math.log(p) + (n - k) * math.log(1 - p)
        if best is None or dev < best[0]:
            best = (dev, tau, a, b)
    return best[1:]


# Sequence context (#35): the called bases two before, one before, at and one
# after a cycle, A, C, G, T only; index = 64 b(-2) + 16 b(-1) + 4 b(0) + b(+1).
BASES = "ACGT"
CONTEXTS = [a + b + c + e for a in BASES for b in BASES for c in BASES for e in BASES]
# Shrinkage of a context factor towards that of the shorter context
# (prev, base, next), and of that towards the base's: a pseudo-count of
# CONTEXT_PRIOR expected calls.
CONTEXT_PRIOR = 20.0


def fit_context(d, out):
    """Per mate, read class and low bin, the factor of each 4-base context.

    From context_classes.tsv (illumina_stats.py context): within the reads
    of a class, the calls of a bin in context x over their expected number
    without a context effect, shrunk towards the factor of the 3-base context
    and of the base, then normalised so that the factors average to 1 over
    the mate's context frequencies (so the per-cycle fractions and class
    levels are unchanged for reads of the sample's composition).
    """
    rows = read_tsv(os.path.join(d, "context_classes.tsv"))
    freq = {}  # (mate, ctx) -> bases
    cell = {}  # (mate, class, ctx) -> [n11, n25, e11, e25]
    for r in rows:
        x = r["prev2"] + r["prev"] + r["base"] + r["next"]
        if any(b not in BASES for b in x):
            continue
        m, k = r["mate"], int(r["class"])
        freq[(m, x)] = freq.get((m, x), 0) + int(r["bases"])
        v = cell.setdefault((m, k, x), [0, 0, 0.0, 0.0])
        v[0] += int(r["n_Q11"])
        v[1] += int(r["n_Q25"])
        v[2] += float(r["exp_Q11"])
        v[3] += float(r["exp_Q25"])
    factors = {}
    for mate in ("1", "2"):
        tot = sum(freq.get((mate, x), 0) for x in CONTEXTS)
        pi = {x: freq.get((mate, x), 0) / tot for x in CONTEXTS}
        for k in range(CLASSES):
            for i, q in enumerate(LOW):
                o = {x: cell.get((mate, k, x), [0, 0, 0, 0])[i] for x in CONTEXTS}
                e = {x: cell.get((mate, k, x), [0, 0, 0, 0])[2 + i] for x in CONTEXTS}
                def pooled(key):
                    po, pe = {}, {}
                    for x in CONTEXTS:
                        po[key(x)] = po.get(key(x), 0) + o[x]
                        pe[key(x)] = pe.get(key(x), 0) + e[x]
                    return po, pe
                o1, e1 = pooled(lambda x: x[2])
                f1 = {b: o1[b] / e1[b] if e1[b] else 1.0 for b in o1}
                o3, e3 = pooled(lambda x: x[1:])
                f3 = {y: (o3[y] + CONTEXT_PRIOR * f1[y[1]]) / (e3[y] + CONTEXT_PRIOR) for y in o3}
                f4 = {x: (o[x] + CONTEXT_PRIOR * f3[x[1:]]) / (e[x] + CONTEXT_PRIOR) for x in CONTEXTS}
                z = sum(pi[x] * f4[x] for x in CONTEXTS)
                for x in CONTEXTS:
                    factors[(mate, k, q, x)] = f4[x] / z
                    out.append((f"context_q{q}", mate, f"{k}:{x}", f4[x] / z))
    tot = sum(freq.values())
    for x in CONTEXTS:
        out.append(("context_freq", "both", x, sum(freq.get((m, x), 0) for m in ("1", "2")) / tot))
    return factors


RAKE_ITERS = 200


def fit_pairs(d, out):
    """The joint class weights of the mates (pair_classes.tsv), raked to the
    per-mate class weights fitted above; None without pair_classes.tsv."""
    path = os.path.join(d, "pair_classes.tsv")
    if not os.path.exists(path):
        return None
    w = {m: [v for p, mate, _, v in out if p == "class_weight" and mate == m] for m in ("1", "2")}
    W = [[0.0] * CLASSES for _ in range(CLASSES)]
    for r in read_tsv(path):
        W[int(r["class_r1"])][int(r["class_r2"])] = float(r["weight"])
    for _ in range(RAKE_ITERS):
        for i in range(CLASSES):
            s = sum(W[i])
            W[i] = [v * w["1"][i] / s for v in W[i]]
        for j in range(CLASSES):
            s = sum(W[i][j] for i in range(CLASSES))
            for i in range(CLASSES):
                W[i][j] *= w["2"][j] / s
    worst = max(abs(sum(W[i][j] for i in range(CLASSES)) - w["2"][j]) for j in range(CLASSES))
    worst = max([worst] + [abs(sum(W[i]) - w["1"][i]) for i in range(CLASSES)])
    print(f"pairs: raked joint class weights, largest margin error {worst:.2g}", file=sys.stderr)
    for i in range(CLASSES):
        for j in range(CLASSES):
            out.append(("class_pair", "pair", f"{i}:{j}", W[i][j]))
    return W


def context_summary(d):
    """Measured quality by called base, by preceding and called base, and
    after GG (context_quals.tsv), per mate and pooled over cycles: bases,
    mean Phred, P(Q11) and P(Q25) with binomial standard errors. Written to
    context_summary.tsv."""
    rows = read_tsv(os.path.join(d, "context_quals.tsv"))
    qcols = [c for c in rows[0] if c.startswith("Q")]
    groups = [
        ("base", lambda r: r["base"]),
        ("prev+base", lambda r: r["prev"] + r["base"]),
        ("base+next", lambda r: r["base"] + r["next"]),
        ("prev2+prev+base", lambda r: r["prev2"] + r["prev"] + r["base"]),
    ]
    table = []
    for mate in ("1", "2"):
        mrows = [r for r in rows if r["mate"] == mate]
        for name, key in groups:
            agg = {}
            for r in mrows:
                x = key(r)
                if any(b not in BASES for b in x):
                    continue
                a = agg.setdefault(x, {c: 0 for c in qcols})
                for c in qcols:
                    a[c] += int(r[c])
            for x in sorted(agg):
                if name == "prev2+prev+base" and x[:2] != "GG":
                    continue
                a = agg[x]
                n = sum(a.values())
                mq = sum(int(c[1:]) * v for c, v in a.items()) / n
                p11, p25 = a["Q11"] / n, a["Q25"] / n
                table.append([mate, name, x, n, f"{mq:.4f}", f"{p11:.6f}", f"{math.sqrt(p11 * (1 - p11) / n):.2g}",
                              f"{p25:.6f}", f"{math.sqrt(p25 * (1 - p25) / n):.2g}"])
    with open(os.path.join(d, "context_summary.tsv"), "w", encoding="ascii") as f:
        f.write("mate\tcontext\tbases_called\tbases\tmean_q\tp_q11\tse_q11\tp_q25\tse_q25\n")
        for r in table:
            f.write("\t".join(str(v) for v in r) + "\n")


def main(d, go_path=None):
    cyc = read_tsv(os.path.join(d, "cycle_quals.tsv"))
    comp = read_tsv(os.path.join(d, "read_composition.tsv"))
    errq = read_tsv(os.path.join(d, "errors_quality.tsv"))
    errc = read_tsv(os.path.join(d, "errors_cycle.tsv"))
    subs = read_tsv(os.path.join(d, "substitutions.tsv"))
    fq = read_tsv(os.path.join(d, "fastq_summary.tsv"))
    trans = read_tsv(os.path.join(d, "transitions.tsv"))
    out = []  # (param, mate, index, value)

    for mate in ("1", "2"):
        rows = [r for r in cyc if r["mate"] == mate]
        L = max(int(r["cycle"]) for r in rows)
        body = [r for r in rows if int(r["cycle"]) < L]
        last = [r for r in rows if int(r["cycle"]) == L][0]
        lam, curves = [], []
        for q in LOW:
            col = f"Q{q}"
            cs = [int(r["cycle"]) for r in body]
            fr = [int(r[col]) / int(r["bases"]) for r in body]
            wt = [int(r["bases"]) for r in body]
            vals = fit_curve(cs, fr, wt)
            for i, v in enumerate(vals):
                out.append((f"low_q{q}", mate, KNOTS[i], v))
            fl = int(last[col]) / int(last["bases"])
            final = fl / interp(L, KNOTS, vals)
            out.append((f"final_factor_q{q}", mate, L, final))
            curves.append([interp(c, KNOTS, vals) for c in range(1, L)] + [final * interp(L, KNOTS, vals)])
            # the mean count per read (the reads' sum of the per-cycle
            # fractions; so for pooled runs of different lengths, #36)
            lam.append(sum(int(r[col]) for r in rows) / max(int(r["bases"]) for r in rows))
        crow = {}
        for r in comp:
            if r["mate"] != mate:
                continue
            k = (int(r["n_Q11"]), int(r["n_Q25"]))
            crow[k] = crow.get(k, 0) + int(r["reads"])
        w, m, ll = em(sorted(crow.items()), lam, CLASSES)
        npar = CLASSES * 3 - 1
        nreads = sum(crow.values())
        print(f"mate {mate}: {CLASSES} classes, log-likelihood {ll:.1f}, BIC {npar * math.log(nreads) - 2 * ll:.1f}", file=sys.stderr)
        out.append(("class_lambda_q11", mate, 0, lam[0]))
        out.append(("class_lambda_q25", mate, 0, lam[1]))
        for k in range(CLASSES):
            out.append(("class_weight", mate, k, w[k]))
            out.append(("class_mult_q11", mate, k, m[k][0]))
            out.append(("class_mult_q25", mate, k, m[k][1]))
        # measured P(b at c+1 | b at c), pooled over cycles
        same, frm = [0, 0], [0, 0]
        for r in trans:
            if r["mate"] != mate:
                continue
            a, b, n = int(r["from"]), int(r["to"]), int(r["count"])
            for i, q in enumerate(LOW):
                if a == q:
                    frm[i] += n
                    same[i] += n * (b == q)
        real = [same[i] / frm[i] for i in range(2)]
        rho = fit_stickiness(w, m, curves, real)
        for i, q in enumerate(LOW):
            out.append((f"stick_q{q}", mate, 0, rho[i]))
            print(f"mate {mate}: P(Q{q}|Q{q}) measured {real[i]:.4f}, model {self_transition(w, m, curves, rho)[i]:.4f}", file=sys.stderr)

    # no-calls: N bases per base, pooled over mates
    ncalls = sum(int(r["value"]) for r in fq if r["key"] == "n_calls")
    bases = sum(int(r["value"]) for r in fq if r["key"] == "bases")
    out.append(("nocall_rate", "both", 0, ncalls / bases))

    # errors
    both = {int(r["phred"]): r for r in errq if r["mate"] == "both"}
    q37 = both[37]
    e37 = int(q37["overlap_errors"]) / int(q37["overlap_bases"])
    variant = int(q37["mismatches"]) / int(q37["aligned"]) - e37
    out.append(("variant_rate", "both", 0, variant))
    out.append(("error_q37", "both", 0, e37))
    for q in LOW:
        agg = {}
        for r in errc:
            if int(r["phred"]) != q:
                continue
            c = int(r["cycle"])
            a = agg.setdefault(c, [0, 0])
            a[0] += int(r["aligned"])
            a[1] += int(r["mismatches"])
        rows = [(c, a, max(k / a - variant, 0)) for c, (a, k) in sorted(agg.items()) if a]
        tau, a, b = fit_error(rows)
        out.append((f"error_q{q}_late", "both", 0, a))
        out.append((f"error_q{q}_early_excess", "both", 0, b))
        out.append((f"error_q{q}_tau", "both", 0, tau))

    # Substitution spectrum, from Q11 mismatches (about 99 % sequencing
    # errors) less the expected sample variants, whose spectrum is that of
    # the Q37 mismatches scaled to the variant share of the Q37 rate.
    by = {(int(r["phred"]), r["ref"]): r for r in subs}
    raw37 = int(q37["mismatches"]) / int(q37["aligned"])
    rate, total = {}, {}
    for x in "ACGT":
        r11, r37 = by[(11, x)], by[(37, x)]
        a11, a37 = int(r11["aligned"]), int(r37["aligned"])
        cnt = {}
        for y in "ACGT":
            if y == x:
                continue
            v = a11 * int(r37["to_" + y]) / a37 * variant / raw37
            cnt[y] = max(int(r11["to_" + y]) - v, 0)
        total[x] = sum(cnt.values())
        rate[x] = total[x] / a11
        for y, c in cnt.items():
            out.append(("sub_matrix", x, y, c / total[x]))
    # Per-base factors, normalised to mean 1 over the human genome's base
    # composition (GRCh38: A, T 29.5 %; C, G 20.5 %).
    comp = {"A": 0.295, "C": 0.205, "G": 0.205, "T": 0.295}
    mean = sum(comp[x] * rate[x] for x in "ACGT")
    for x in "ACGT":
        out.append(("base_factor", x, 0, rate[x] / mean))

    fit_pairs(d, out)

    have_context = os.path.exists(os.path.join(d, "context_classes.tsv"))
    if have_context:
        fit_context(d, out)
        context_summary(d)
    elif go_path:
        raise SystemExit("context_classes.tsv missing: run illumina_stats.py context first")

    with open(os.path.join(d, "fit.tsv"), "w", encoding="ascii") as f:
        f.write("param\tmate\tindex\tvalue\n")
        for p, mate, i, v in out:
            f.write(f"{p}\t{mate}\t{i}\t{v:.6g}\n")
    if go_path:
        info = read_profile(d)
        with open(go_path, "w", encoding="ascii") as f:
            f.write(go_source(out, info))


def read_profile(d):
    """The run description of RESULTS_DIR/profile.tsv (written by
    illumina_calibrate.sh) and the headline statistics of the measured reads
    (fastq_summary.tsv, sam_summary.tsv)."""
    path = os.path.join(d, "profile.tsv")
    if not os.path.exists(path):
        raise SystemExit(f"{path} missing: illumina_calibrate.sh writes it (name, accession, ...)")
    info = {r["key"]: r["value"] for r in read_tsv(path)}
    for k in ("name", "accession", "instrument", "center", "library", "study", "pairs"):
        if k not in info:
            raise SystemExit(f"{path}: no {k}")
    # A pooled profile (several runs, illumina_calibrate.sh): the runs' ENA
    # metadata, one row per run.
    info["runs"] = []
    if "id" in info:
        info["runs"] = read_tsv(os.path.join(d, "profile_runs.tsv"))
    fq = {(r["mate"], r["key"]): r["value"] for r in read_tsv(os.path.join(d, "fastq_summary.tsv"))}
    sam = {r["key"]: r["value"] for r in read_tsv(os.path.join(d, "sam_summary.tsv"))}
    info["stats"] = {
        "MeanQ": [float(fq[(m, "mean_q")]) for m in "12"],
        "FracQ30": [float(fq[(m, "frac_ge_q30")]) for m in "12"],
        "FracQ11": [float(fq.get((m, "frac_Q11"), 0)) for m in "12"],
        "FracQ25": [float(fq.get((m, "frac_Q25"), 0)) for m in "12"],
        "ReadLength": int(fq[("1", "max_length")]),
        "InsertMedian": float(sam["insert_median"]),
        "InsertMean": float(sam.get("insert_core_mean", sam["insert_mean"])),
        "InsertSD": float(sam.get("insert_core_sd", sam["insert_sd"])),
    }
    return info


def go_quote(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def round5(x):
    """x rounded to the nearest multiple of 5 (the insert-size defaults)."""
    return int(5 * math.floor(x / 5 + 0.5))


def go_source(out, info):
    """Render the fitted parameters of one run as a Go file holding its
    Profile (internal/genomics/illumina/profile_<accession>.go)."""
    P = {}
    for p, mate, i, v in out:
        P[(p, str(mate), str(i))] = float(f"{v:.6g}")

    def g(p, mate="both", i="0"):
        return repr(P[(p, str(mate), str(i))])

    acc = info["accession"]
    st = info["stats"]
    L = st["ReadLength"]
    sfx = info.get("id", acc).upper()
    runs = info["runs"]
    if runs:
        about = [
            f"// profile{sfx} is fitted to {int(info['pairs']):,} read pairs pooled from {len(runs)} {info['instrument']}",
            f"// runs (human WGS), the first {int(info['pairs_per_run']):,} pairs of each of ENA runs",
            f"// {', '.join(r['accession'] for r in runs)} (#36); see",
        ]
    else:
        about = [
            f"// profile{sfx} is fitted to the first {int(info['pairs']):,} read pairs (2x{L}) of ENA",
            f"// run {acc} ({info['instrument']}, human WGS); see",
        ]
    lines = [
        "// Code generated by scripts/acceptance/calibrate/illumina/fit_illumina.py; DO NOT EDIT.",
        "",
        "package illumina",
        "",
        'import "github.com/scttfrdmn/phantom-science/pkg/types"',
        "",
    ] + about + [
        "// scripts/acceptance/illumina_calibrate.sh and the measured tables in",
        f"// {info.get('results', 'scripts/acceptance/calibrate/illumina/results')}.",
        f"var profile{sfx} = Profile{{",
        f"	Name: {go_quote(info['name'])},",
        "	Run: RunInfo{",
        f"		Accession: {go_quote(acc)}, Study: {go_quote(info['study'])},",
        f"		Instrument: {go_quote(info['instrument'])},",
        f"		Center: {go_quote(info['center'])},",
        f"		Library: {go_quote(info['library'])},",
        f"		ReadLength: {L}, Pairs: {int(info['pairs'])},",
        "	},",
    ]
    if runs:
        lines.append("	Runs: []RunInfo{")
        for r in runs:
            lines += [
                "		{",
                f"			Accession: {go_quote(r['accession'])}, Study: {go_quote(r['study'])},",
                f"			Instrument: {go_quote(r['instrument'])},",
                f"			Center: {go_quote(r['center'])},",
                f"			Library: {go_quote(r['library'])},",
                f"			ReadLength: {max(int(x) for x in r['read_length'].split(','))}, Pairs: {int(r['pairs'])},",
                "		},",
            ]
        lines.append("	},")
    lines += [
        "	Stats: RunStats{",
        f"		MeanQ: [2]float64{{{st['MeanQ'][0]:.4g}, {st['MeanQ'][1]:.4g}}}, FracQ30: [2]float64{{{st['FracQ30'][0]:.4g}, {st['FracQ30'][1]:.4g}}},",
        f"		FracQ11: [2]float64{{{st['FracQ11'][0]:.4g}, {st['FracQ11'][1]:.4g}}}, FracQ25: [2]float64{{{st['FracQ25'][0]:.4g}, {st['FracQ25'][1]:.4g}}},",
        f"		InsertMedian: {st['InsertMedian']:.6g}, InsertMean: {st['InsertMean']:.6g}, InsertSD: {st['InsertSD']:.6g},",
        "	},",
    ]
    for mate in ("1", "2"):
        final = [i for (p, m, i) in P if p == "final_factor_q11" and m == mate][0]
        lines.append(f"	Read{mate}: QualityParams{{")
        lines.append("		Knots: []int{" + ", ".join(str(k) for k in KNOTS) + "},")
        lines.append("		Low:   []float64{" + ", ".join(g("low_q11", mate, k) for k in KNOTS) + "},")
        lines.append("		Mid:   []float64{" + ", ".join(g("low_q25", mate, k) for k in KNOTS) + "},")
        lines.append(f"		FinalLow: {g('final_factor_q11', mate, final)}, FinalMid: {g('final_factor_q25', mate, final)},")
        lines.append("		Classes: []ReadClass{")
        for k in range(CLASSES):
            lines.append(f"			{{Weight: {g('class_weight', mate, k)}, Low: {g('class_mult_q11', mate, k)}, Mid: {g('class_mult_q25', mate, k)}, Context: read{mate}Context{sfx}[{k}]}},")
        lines.append("		},")
        lines.append(f"		StickLow: {g('stick_q11', mate)}, StickMid: {g('stick_q25', mate)},")
        lines.append(f"		NoCall: {g('nocall_rate')},")
        lines.append("	},")
    if ("class_pair", "pair", "0:0") in P:
        lines.append("	Pairs: PairParams{Joint: [][]float64{")
        for i in range(CLASSES):
            lines.append("		{" + ", ".join(g("class_pair", "pair", f"{i}:{j}") for j in range(CLASSES)) + "},")
        lines.append("	}},")
    early = lambda q: repr(float(f"{P[(f'error_q{q}_late', 'both', '0')] + P[(f'error_q{q}_early_excess', 'both', '0')]:.6g}"))
    lines += [
        "	Errors: ErrorParams{",
        f"		High:     {g('error_q37')},",
        f"		MidEarly: {early(25)}, MidLate: {g('error_q25_late')}, MidTau: {g('error_q25_tau')},",
        f"		LowEarly: {early(11)}, LowLate: {g('error_q11_late')}, LowTau: {g('error_q11_tau')},",
        "		BaseFactor: [4]float64{" + ", ".join(g("base_factor", x) for x in "ACGT") + "},",
        "		Matrix: [4][4]float64{",
    ]
    for x in "ACGT":
        lines.append("			{" + ", ".join("0" if x == y else g("sub_matrix", x, y) for y in "ACGT") + f"}}, // {x}")
    lines += [
        "		},",
        "	},",
        "	// The measured insert sizes of proper pairs, lognormal, mean and SD",
        "	// rounded to 5.",
        f"	InsertSize: types.InsertSizeSpec{{Dist: types.InsertSizeLognormal, Mean: {round5(st['InsertMean'])}, SD: {round5(st['InsertSD'])}}},",
        f"	freq: &contextFreq{sfx},",
        "}",
        "",
    ]
    lines += context_source(out, sfx)
    return "\n".join(lines)


def context_source(out, sfx):
    """The context frequencies and factors (fit_context) as Go source."""
    C = {}
    for p, mate, i, v in out:
        if p.startswith("context_"):
            C[(p, str(mate), str(i))] = v

    def row(vals):
        # 8 values per line
        return ["		" + ", ".join(f"{v:.4g}" for v in vals[j : j + 8]) + "," for j in range(0, len(vals), 8)]

    lines = [
        f"// contextFreq{sfx} is the frequency of each 4-base context (index 64 b(-2) +",
        "// 16 b(-1) + 4 b(0) + b(+1), A, C, G, T = 0-3) in the measured reads, both",
        "// mates; the factors average to 1 over it.",
        f"var contextFreq{sfx} = [256]float64{{",
    ]
    lines += row([C[("context_freq", "both", x)] for x in CONTEXTS])
    lines += ["}", ""]
    for mate in ("1", "2"):
        lines.append(f"// read{mate}Context{sfx}[k] are the context factors of read {mate}'s class k: P(Q11)")
        lines.append("// and P(Q25) at a cycle are multiplied by Low and Mid of the cycle's context.")
        lines.append(f"var read{mate}Context{sfx} = [...]*ContextFactors{{")
        for k in range(CLASSES):
            lines.append(f"	NewContextFactors(&contextFreq{sfx},")
            for q in LOW:
                lines.append("		[256]float64{")
                lines += ["	" + l for l in row([C[(f"context_q{q}", mate, f"{k}:{x}")] for x in CONTEXTS])]
                lines.append("		}," if q == LOW[0] else "		}),")
        lines += ["}", ""]
    return lines


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
