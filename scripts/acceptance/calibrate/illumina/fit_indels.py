"""Fit the indel-error model of internal/genomics/illumina to the tables of
`illumina_stats.py indels` (#30).

Run with `python3 -I fit_indels.py RESULTS GO_FILE SHAPE NAME=DIR... [@DIR...] [%DIR=NAME]`:
RESULTS is the directory of the measured runs (one subdirectory each),
SHAPE the subdirectory the shape of the model is fitted to (the pooled
runs), each NAME=DIR a profile and the subdirectory of the run (or pool)
it was fitted to, and each @DIR a held-out run, compared with the shape's
expectation (not fitted). Writes RESULTS/indel_fit.tsv and GO_FILE (the
profiles' IndelParams, gofmt'd by the caller).

Model (see IndelParams in internal/genomics/illumina): an indel error is
a slip in a homopolymer run (per run, by its length h), a slip of whole
units in a short tandem repeat (units of 2-6 bases, per tract, by its
copies c), or a background insertion or deletion (per base). For each kind
(insertion, deletion) the rate per run is exp(a + b (h - 1)) and per tract
exp(a + b (c - 2)), fitted by Poisson maximum likelihood to the events and
opportunities of SHAPE (indel_classes.tsv), capped at RATE_CAP and constant
past the last measured length; the background rate is the events per base.
The fits keep each class's total: the expected events of SHAPE equal the
observed. The sizes (bases of a homopolymer slip, apart for runs shorter
and longer than LONG_RUN; units of a tandem-repeat slip; bases of a
background indel) and the qualities of inserted bases and of the bases
either side of a deletion are SHAPE's. A profile scales the insertion and
the deletion rates by its run's observed events over those the shape
expects on the run's opportunities, so it reproduces its run's rates.

The measurement does not find every error (errors in long repeats and
long indels are lost to the aligner more often). With %DIR=PROFILE, DIR
holds reads generated with PROFILE's model fitted without this step and
measured the same way (the closed loop): the insertion and the deletion
rates are multiplied by PROFILE's run's measured rate per base over the
generated reads' (indel_fit.tsv, "recovery", with the errors generated and
found by group for the record: short homopolymer runs, longer repeats,
background), so that reads generated with the final model measure like
the real ones.

indel_fit.tsv also tests whether SHAPE's rates depend on the cycle or the
mate ("dependence"); the model has neither.
"""

import math
import os
import sys

HP_CAP, STR_CAP = 30, 20  # illumina_stats.py
LONG_RUN = 10  # homopolymer runs of this length or more have their own sizes
SIZE_CAP = 30  # sizes are tabulated up to this many bases (units)
RATE_CAP = 0.25
RECOV_CLAMP = 2.0
BINS = (2, 11, 25, 37)


def read_tsv(path):
    with open(path, encoding="ascii") as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f]


def poisson_ci(k, conf=0.95):
    a = (1 - conf) / 2

    def cdf(k, mu):
        if k < 0:
            return 0.0
        if mu <= 0:
            return 1.0
        lm = math.log(mu)
        return min(1.0, sum(math.exp(i * lm - mu - math.lgamma(i + 1)) for i in range(k + 1)))

    def solve(f, lo, hi):
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid):
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    hi0 = k + 20 * math.sqrt(k + 1) + 20
    lower = 0.0 if k == 0 else solve(lambda mu: 1 - cdf(k - 1, mu) < a, 0.0, hi0)
    return lower, solve(lambda mu: cdf(k, mu) > a, 0.0, hi0)


def classes(d):
    """{(op, class, x): (opportunities, events)} of a run's indel_classes.tsv
    (x the run length or copies, 0 for background)."""
    out = {}
    for r in read_tsv(os.path.join(d, "indel_classes.tsv")):
        x = 0 if r["class"] == "bg" else int(r["bin"].rstrip("+"))
        out[(r["op"], r["class"], x)] = (int(r["opportunities"]), int(r["events"]))
    return out


def fit_loglinear(points, x0):
    """Poisson MLE of log rate = a + b (x - x0) for points [(x, n, k)]
    (damped Newton steps)."""
    pts = [(x - x0, n, k) for x, n, k in points if n > 0]

    def loglik(a, b):
        return sum(k * (a + b * z) - n * math.exp(min(a + b * z, 50)) for z, n, k in pts)

    a = math.log(max(sum(k for _, _, k in pts), 0.5) / sum(n for _, n, _ in pts))
    b = 0.0
    for _ in range(200):
        g = [0.0, 0.0]
        h = [[0.0, 0.0], [0.0, 0.0]]
        for z, n, k in pts:
            mu = n * math.exp(min(a + b * z, 50))
            g[0] += k - mu
            g[1] += (k - mu) * z
            h[0][0] += mu
            h[0][1] += mu * z
            h[1][1] += mu * z * z
        det = h[0][0] * h[1][1] - h[0][1] ** 2
        if det <= 0:
            break
        da = (h[1][1] * g[0] - h[0][1] * g[1]) / det
        db = (h[0][0] * g[1] - h[0][1] * g[0]) / det
        t, base = 1.0, loglik(a, b)
        while t > 1e-6 and loglik(a + t * da, b + t * db) < base:
            t /= 2
        a += t * da
        b += t * db
        if abs(t * da) < 1e-12 and abs(t * db) < 1e-12:
            break
    return a, b


def sizes(d, op, cls, sel):
    """The size distribution (P(size = 1..SIZE_CAP)) of a class's events."""
    n = [0] * SIZE_CAP
    for r in read_tsv(os.path.join(d, "indel_events.tsv")):
        if r["op"] == op and r["class"] == cls and sel(int(r["copies"])):
            n[min(int(r["size"]), SIZE_CAP) - 1] += 1
    tot = sum(n)
    while n and n[-1] == 0:
        n.pop()
    return [x / tot for x in n] if tot else [1.0], tot


def quals(d):
    rows = [r for r in read_tsv(os.path.join(d, "indel_quality.tsv")) if r["mate"] == "both"]
    out = {}
    for col in ("inserted_bases", "deletion_flank_bases"):
        n = {int(r["phred"]): int(r[col]) for r in rows}
        tot = sum(n.values())
        out[col] = [n.get(q, 0) / tot for q in BINS]
    return out


def shape(d):
    """The fitted shape: rates per (op, class, x) and the fit rows."""
    c = classes(d)
    model, rows = {}, []
    for op in "ID":
        for cls, cap, x0 in (("hp", HP_CAP, 1), ("str", STR_CAP, 2)):
            pts = [(x, n, k) for (o, cl, x), (n, k) in c.items() if o == op and cl == cls]
            a, b = fit_loglinear(pts, x0)
            for x in range(x0, cap + 1):
                model[(op, cls, x)] = min(RATE_CAP, math.exp(a + b * (x - x0)))
            obs = sum(k for _, _, k in pts)
            exp_ = sum(n * model[(op, cls, x)] for x, n, _ in pts)
            dev = 2 * sum((k * math.log(k / (n * model[(op, cls, x)])) if k else 0) - (k - n * model[(op, cls, x)]) for x, n, k in pts if n)
            rows += [[op, cls, "intercept_a", f"{a:.6g}"], [op, cls, "slope_b", f"{b:.6g}"],
                     [op, cls, "observed", obs], [op, cls, "expected", f"{exp_:.6g}"],
                     [op, cls, "deviance", f"{dev:.4g}"], [op, cls, "points", sum(1 for _, n, _ in pts if n)]]
        n, k = c[(op, "bg", 0)]
        model[(op, "bg", 0)] = k / n
        rows += [[op, "bg", "rate", f"{k / n:.6g}"], [op, "bg", "observed", k]]
    return model, rows


def expected(model, c, op):
    return sum(n * model[(o, cl, x)] for (o, cl, x), (n, _) in c.items() if o == op)


def chi2_sf_even(x, df):
    """P(X > x) for a chi-square variable with an even df."""
    t, term, tot = x / 2, 1.0, 1.0
    for i in range(1, df // 2):
        term *= t / i
        tot += term
    return math.exp(-t) * tot


def dependence(d):
    """Rows testing whether the rates depend on the cycle (Pearson's
    chi-square of the cycle blocks against one rate, 4 degrees of freedom)
    and on the mate (a binomial split of the events, the mates reading the
    same overlap)."""
    rows = []
    blocks = [r for r in read_tsv(os.path.join(d, "indel_cycles.tsv")) if r["mate"] == "both"]
    for op, col in (("I", "ins"), ("D", "del")):
        n = [int(r["bases"]) for r in blocks]
        k = [int(r[col]) for r in blocks]
        rate = sum(k) / sum(n)
        chi2 = sum((ki - ni * rate) ** 2 / (ni * rate) for ki, ni in zip(k, n) if ni)
        df = len(blocks) - 1
        p = chi2_sf_even(chi2, df) if df % 2 == 0 else float("nan")
        rows.append(["dependence", op, "cycle", "chi2", f"chi2={chi2:.3g} df={df} p={p:.3g} blocks={','.join(r['cycles'] for r in blocks)} "
                     f"events={','.join(map(str, k))}"])
    summ = {r["key"]: r for r in read_tsv(os.path.join(d, "indel_summary.tsv"))}
    for op, lab in (("I", "ins"), ("D", "del")):
        k1, k2 = int(summ[f"mate1_{lab}_rate_per_base"]["events"]), int(summ[f"mate2_{lab}_rate_per_base"]["events"])
        n = k1 + k2
        # two-sided exact binomial test of k1 ~ Bin(n, 1/2)
        pk = [math.comb(n, i) / 2 ** n for i in range(n + 1)]
        p = min(1.0, sum(x for x in pk if x <= pk[k1] * (1 + 1e-9)))
        rows.append(["dependence", op, "mate", "binomial", f"mate1={k1} mate2={k2} p={p:.3g}"])
    return rows


def group(cls, x):
    """The recovery group of a class: short homopolymer runs, longer
    repeats (runs of LONG_RUN or more, tandem repeats), background."""
    if cls == "hp" and x < LONG_RUN:
        return "short"
    return "bg" if cls == "bg" else "repeat"


GROUPS = ("short", "repeat", "bg")


def go_floats(xs, typed=True):
    return ("[]float64" if typed else "") + "{" + ", ".join(f"{x:.6g}" for x in xs) + "}"


def main(argv):
    if len(argv) < 5:
        raise SystemExit(__doc__)
    results, go_path, shape_dir = argv[1], argv[2], argv[3]
    profs = [a.split("=", 1) for a in argv[4:] if a[0] not in "@%"]
    checks = [a[1:] for a in argv[4:] if a.startswith("@")]
    recov = [a[1:].split("=", 1) for a in argv[4:] if a.startswith("%")]
    sd = os.path.join(results, shape_dir)
    model, rows = shape(sd)
    out = [["shape", shape_dir, "", "", ""]]
    out += [["fit", r[0], r[1], r[2], r[3]] for r in rows]
    sz = {}
    for op in "ID":
        sz[(op, "hp", 0)], n0 = sizes(sd, op, "hp", lambda c: c < LONG_RUN)
        sz[(op, "hp", 1)], n1 = sizes(sd, op, "hp", lambda c: c >= LONG_RUN)
        sz[(op, "str", 0)], n2 = sizes(sd, op, "str", lambda c: True)
        sz[(op, "bg", 0)], n3 = sizes(sd, op, "bg", lambda c: True)
        out += [["sizes_events", op, cls, "", n] for cls, n in (("hp<10", n0), ("hp10+", n1), ("str", n2), ("bg", n3))]
    out += dependence(sd)
    q = quals(sd)
    out += [["quality", "inserted", BINS[i], "", f"{v:.6g}"] for i, v in enumerate(q["inserted_bases"])]
    out += [["quality", "deletion_flank", BINS[i], "", f"{v:.6g}"] for i, v in enumerate(q["deletion_flank_bases"])]

    # Runs no profile is fitted to: their measured errors against the
    # shape's (the pooled profile's) expectation on their opportunities.
    for sub in checks:
        c = classes(os.path.join(results, sub))
        for op in "ID":
            obs = sum(k for (o, _, _), (_, k) in c.items() if o == op)
            exp_ = expected(model, c, op)
            lo, hi = poisson_ci(obs)
            bases = c[(op, "bg", 0)][0]
            out.append(["heldout", sub, "", op, f"observed={obs} expected={exp_:.6g} ratio={obs / exp_:.4g} "
                        f"ci=[{lo:.6g},{hi:.6g}] expected_in_ci={'yes' if lo <= exp_ <= hi else 'no'} rate={obs / bases:.6g}"])

    def scales(sub):
        c = classes(os.path.join(results, sub))
        out = {}
        for op in "ID":
            obs = sum(k for (o, _, _), (_, k) in c.items() if o == op)
            exp_ = expected(model, c, op)
            out[op] = obs / exp_ if exp_ else 1.0
        return out

    # Recovery (the closed loop): reads generated with a profile's model as
    # fitted so far, measured the same way, show what share of each group's
    # errors the measurement finds. The rates are divided by it, so that
    # generated reads measure like the real ones. RECOV_CLAMP bounds it.
    factor = {(op, g): 1.0 for op in "ID" for g in GROUPS}
    subs = dict(profs)
    for rsub, rprof in recov:
        sc = scales(subs[rprof])
        c = classes(os.path.join(results, rsub))
        for op in "ID":
            # By group for the record; the correction is one factor per kind
            # (insertions, deletions): errors found in another group than
            # they were made in (a background deletion measured as a
            # homopolymer slip) make per-group factors overshoot.
            for g in GROUPS:
                obs = sum(k for (o, cl, x), (n, k) in c.items() if o == op and group(cl, x) == g)
                exp_ = sum(n * model[(o, cl, x)] * sc[op] for (o, cl, x), (n, k) in c.items() if o == op and group(cl, x) == g)
                out.append(["recovery", op, g, "detail", f"set={rsub} profile={rprof} found={obs} generated={exp_:.6g}"])
            obs = sum(k for (o, _, _), (_, k) in c.items() if o == op)
            exp_ = sum(n * model[(o, cl, x)] * sc[op] for (o, cl, x), (n, _) in c.items() if o == op)
            gen_rate = obs / c[(op, "bg", 0)][0]
            rc = classes(os.path.join(results, subs[rprof]))
            real = sum(k for (o, _, _), (_, k) in rc.items() if o == op)
            real_rate = real / rc[(op, "bg", 0)][0]
            # The factor: the profile's run's measured rate per base over
            # that of its generated reads, so that generated reads measure
            # like the run (it also covers what the measurement loses).
            f = min(max(real_rate / gen_rate, 1 / RECOV_CLAMP), RECOV_CLAMP) if obs else 1.0
            out.append(["recovery", op, "all", "detail", f"set={rsub} profile={rprof} found={obs} generated={exp_:.6g} "
                        f"rate={gen_rate:.6g} run_rate={real_rate:.6g}"])
            for g in GROUPS:
                factor[(op, g)] = f
                out.append(["recovery", op, g, "factor", f"{f:.6g}"])
    measured = model
    model = {k: min(RATE_CAP, v * factor[(k[0], group(k[1], k[2]))]) for k, v in measured.items()}

    go = []
    for name, sub in profs:
        c = classes(os.path.join(results, sub))
        scale = {}
        for op in "ID":
            obs = sum(k for (o, _, _), (_, k) in c.items() if o == op)
            exp_ = expected(measured, c, op)
            scale[op] = obs / exp_ if exp_ else 1.0
            bases = c[(op, "bg", 0)][0]
            lo, hi = poisson_ci(obs)
            out.append(["profile", name, sub, op, f"observed={obs} expected={exp_:.6g} scale={scale[op]:.6g} "
                        f"rate={obs / bases:.6g} ci=[{lo / bases:.6g},{hi / bases:.6g}] bases={bases}"])

        def table(op, cls, x0, cap):
            return go_floats([float(f"{model[(op, cls, x)] * scale[op]:.6g}") for x in range(x0, cap + 1)])

        var = "profilePOOLED" if name.endswith("-pooled") else "profile" + sub.upper()
        go.append(f"""	// {name}: {sub}, insertions ×{scale['I']:.4g} and deletions ×{scale['D']:.4g} of the shape.
	{var}.Indels = IndelParams{{
		HPIns:  {table('I', 'hp', 1, HP_CAP)},
		HPDel:  {table('D', 'hp', 1, HP_CAP)},
		STRIns: {table('I', 'str', 2, STR_CAP)},
		STRDel: {table('D', 'str', 2, STR_CAP)},
		BgIns:  {model[('I', 'bg', 0)] * scale['I']:.6g},
		BgDel:  {model[('D', 'bg', 0)] * scale['D']:.6g},
		HPInsSize:  [2][]float64{{{go_floats(sz[('I', 'hp', 0)], False)}, {go_floats(sz[('I', 'hp', 1)], False)}}},
		HPDelSize:  [2][]float64{{{go_floats(sz[('D', 'hp', 0)], False)}, {go_floats(sz[('D', 'hp', 1)], False)}}},
		STRInsSize: {go_floats(sz[('I', 'str', 0)])},
		STRDelSize: {go_floats(sz[('D', 'str', 0)])},
		BgInsSize:  {go_floats(sz[('I', 'bg', 0)])},
		BgDelSize:  {go_floats(sz[('D', 'bg', 0)])},
		InsQual: [4]float64{{{", ".join(f"{v:.6g}" for v in q["inserted_bases"])}}},
		DelQual: [4]float64{{{", ".join(f"{v:.6g}" for v in q["deletion_flank_bases"])}}},
	}}""")
    with open(os.path.join(results, "indel_fit.tsv"), "w", encoding="ascii") as f:
        f.write("section\tkey\tclass\tparameter\tvalue\n")
        for r in out:
            f.write("\t".join(str(x) for x in r) + "\n")
    recov_note = ""
    if recov:
        recov_note = (", and multiplied by the run's measured rate over that of reads\n"
                      f"// generated with the model (insertions ×{factor[('I', 'bg')]:.4g}, deletions ×{factor[('D', 'bg')]:.4g})")
    with open(go_path, "w", encoding="utf-8") as f:
        f.write(f"""// Code generated by scripts/acceptance/calibrate/illumina/fit_indels.py from
// scripts/acceptance/calibrate/illumina/results/indels; DO NOT EDIT.

package illumina

// The indel-error models of the profiles (#30): the shape fitted to
// {shape_dir} (the fitted runs pooled), each profile's insertion and deletion rates
// scaled to its run's measured ones{recov_note}
// (results/indels/indel_fit.tsv).
func init() {{
""" + "\n".join(go) + "\n}\n")


if __name__ == "__main__":
    main(sys.argv)
