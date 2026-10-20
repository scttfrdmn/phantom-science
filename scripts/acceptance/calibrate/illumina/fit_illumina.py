"""Fit the phantom-science Illumina model to the tables of illumina_stats.py (#27).

Run with `python3 -I fit_illumina.py RESULTS_DIR [GO_FILE]`; writes
RESULTS_DIR/fit.tsv (parameter, mate, index, value) and, with GO_FILE, the
parameters as Go source (internal/genomics/illumina/calibration.go).

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
* Substitution rate per bin and cycle: the mismatch rate against hg38 minus
  the rate of sample variants. The variant rate is the Q37 mismatch rate
  minus the Q37 error rate measured where mates overlap (a mismatch in one
  mate only). Q37 has one rate (its mismatch rate does not change with
  cycle); Q11 and Q25 follow e(c) = A + B exp(-(c-1)/tau), fitted by
  weighted least squares with tau on a grid.
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


def fit_error(rows_by_cycle):
    """Fit e(c) = A + B exp(-(c-1)/tau) to (cycle, aligned, rate) rows."""
    best = None
    for tau in [2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300]:
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
            lam.append(sum(int(r[col]) / int(r["bases"]) for r in rows))
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

    with open(os.path.join(d, "fit.tsv"), "w", encoding="ascii") as f:
        f.write("param\tmate\tindex\tvalue\n")
        for p, mate, i, v in out:
            f.write(f"{p}\t{mate}\t{i}\t{v:.6g}\n")
    if go_path:
        with open(go_path, "w", encoding="ascii") as f:
            f.write(go_source(out))


def go_source(out):
    """Render the fitted parameters as internal/genomics/illumina/calibration.go."""
    P = {}
    for p, mate, i, v in out:
        P[(p, str(mate), str(i))] = float(f"{v:.6g}")

    def g(p, mate="both", i="0"):
        return repr(P[(p, str(mate), str(i))])

    lines = [
        "// Code generated by scripts/acceptance/calibrate/illumina/fit_illumina.py; DO NOT EDIT.",
        "",
        "package illumina",
        "",
        "// Parameters fitted to the first 200,000 read pairs (2x151) of ENA run",
        "// DRR276626 (Illumina NovaSeq 6000, human WGS); see",
        "// scripts/acceptance/illumina_calibrate.sh and the measured tables in",
        "// scripts/acceptance/calibrate/illumina/results.",
        "",
        "// novaSeqKnots are the cycles (1-based) of the piecewise-linear low-bin",
        "// curves.",
        "var novaSeqKnots = []int{" + ", ".join(str(k) for k in KNOTS) + "}",
        "",
    ]
    for mate in ("1", "2"):
        lines.append(f"// Read{mate} is the quality model of read {mate} (of single-end reads too)." if mate == "1" else "// Read2 is the quality model of read 2.")
        lines.append(f"var Read{mate} = QualityParams{{")
        lines.append("	Knots: novaSeqKnots,")
        lines.append("	Low:   []float64{" + ", ".join(g("low_q11", mate, k) for k in KNOTS) + "},")
        lines.append("	Mid:   []float64{" + ", ".join(g("low_q25", mate, k) for k in KNOTS) + "},")
        lines.append(f"	FinalLow: {g('final_factor_q11', mate, 151)}, FinalMid: {g('final_factor_q25', mate, 151)},")
        lines.append("	Classes: []ReadClass{")
        for k in range(CLASSES):
            lines.append(f"		{{Weight: {g('class_weight', mate, k)}, Low: {g('class_mult_q11', mate, k)}, Mid: {g('class_mult_q25', mate, k)}}},")
        lines.append("	},")
        lines.append(f"	StickLow: {g('stick_q11', mate)}, StickMid: {g('stick_q25', mate)},")
        lines.append(f"	NoCall: {g('nocall_rate')},")
        lines.append("}")
        lines.append("")
    early = lambda q: repr(float(f"{P[(f'error_q{q}_late', 'both', '0')] + P[(f'error_q{q}_early_excess', 'both', '0')]:.6g}"))
    lines += [
        "// NovaSeqErrors is the substitution model measured against GRCh38.",
        "var NovaSeqErrors = ErrorParams{",
        f"	High:     {g('error_q37')},",
        f"	MidEarly: {early(25)}, MidLate: {g('error_q25_late')}, MidTau: {g('error_q25_tau')},",
        f"	LowEarly: {early(11)}, LowLate: {g('error_q11_late')}, LowTau: {g('error_q11_tau')},",
        "	BaseFactor: [4]float64{" + ", ".join(g("base_factor", x) for x in "ACGT") + "},",
        "	Matrix: [4][4]float64{",
    ]
    for x in "ACGT":
        lines.append("		{" + ", ".join("0" if x == y else g("sub_matrix", x, y) for y in "ACGT") + f"}}, // {x}")
    lines += ["	},", "}", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
