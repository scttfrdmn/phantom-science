"""Held-out validation of the Illumina model (#35).

Run with `python3 -I` (no third-party packages). Subcommands:

  basequal OUTDIR R1.fastq R2.fastq
      Per mate and called base: bases, mean quality and the fraction of
      each quality value (base_quality.tsv); how the two mates' qualities
      go together: the correlation of the mates' low-quality (Phred < 20,
      and Phred < 30) base counts, the pairs with both mates poor against
      the number expected if the mates were independent, and the mean
      correlation of the mates' 10 cycle-block mean qualities (the pair
      classifier's feature) (pair_coupling.tsv).

  classifiers OUTDIR STATS_ROOT PROFILE=RUN,... -- PROFILE@RUN=E2E.json ... -- A:B=E2E.json ...
      Reads run_e2e_tests reports: one per profile and run (real vs real,
      and real vs the profile's generated reads), and one per pair of runs
      (-real-fastq-b, real run A vs real run B). PROFILE=RUN names the run
      each profile was fitted to; the other runs are held out. Writes
      classifier_matrix.tsv (every comparison), classifier_spread.tsv (per
      read and feature set, min / median / max of: the real-vs-real nulls,
      a profile against its own run, against the other fitted runs, against
      the held-out runs, and real run against real run), profile_matrix.tsv
      (the accuracy of every profile's reads and every real run against
      every run) and profile_matching.tsv: per run, the profile chosen by
      the matching rule (the profile whose run is closest in the mean
      quality of read 1 and read 2, STATS_ROOT/RUN/fastq_summary.tsv; for a
      fitted run, among the other profiles), its accuracy, the best other
      profile's, and the accuracy of the other real runs against the run.

  calibration OUTDIR STATS_ROOT PROFILE=RUN,... SET ...
      Compares calibration statistics of each read set (illumina_stats.py
      fastq and sam output, and basequal output, in STATS_ROOT/SET; the
      generated reads of a profile are the set generated-PROFILE). Writes
      calibration_stats.tsv (one row per statistic and mate, one column per
      read set, the spread of the held-out runs, and each profile's
      generated reads minus its run) and cycle_delta.tsv (mean |difference|
      of the per-cycle quality fractions between every two read sets).

All output is tab-separated with a header line.
"""

import json
import math
import os
import statistics
import sys

LOW_Q = 20
POOR_PAIR = 5  # a mate is poor with at least this many bases below LOW_Q
LOW_Q30 = 30
BLOCKS = 10  # cycle blocks of the block-mean correlation (internal/testing.PairFeatures)


def block_means(q):
    n = len(q)
    s, c = [0] * BLOCKS, [0] * BLOCKS
    for i, ch in enumerate(q):
        b = i * BLOCKS // n
        s[b] += ord(ch) - 33
        c[b] += 1
    return [x / k if k else 0.0 for x, k in zip(s, c)]


def pearson(x, y):
    """As internal/testing.pearson: 0 when either is constant."""
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx < 1e-12 or syy < 1e-12:
        return 0.0
    return sxy / math.sqrt(sxx * syy)


def write_tsv(path, header, rows):
    with open(path, "w", encoding="ascii") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def read_tsv(path):
    with open(path, encoding="ascii") as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f if line.strip()]


def fmt(x, digits=6):
    if x is None:
        return ""
    if isinstance(x, int):
        return str(x)
    return f"{x:.{digits}g}"


def fastq_records(path):
    with open(path, "r", encoding="ascii") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().rstrip("\n")
            f.readline()
            q = f.readline().rstrip("\n")
            yield s, q


# ---------------------------------------------------------------- basequal


def basequal(outdir, r1, r2):
    rows = []
    lows = [[], []]
    lows30 = [[], []]
    blocks = [[], []]
    for mate, path in enumerate((r1, r2), start=1):
        counts = {}  # base -> {qual char: n}
        for s, q in fastq_records(path):
            low = low30 = 0
            for b, ch in zip(s, q):
                d = counts.setdefault(b, {})
                d[ch] = d.get(ch, 0) + 1
                if ord(ch) - 33 < LOW_Q:
                    low += 1
                if ord(ch) - 33 < LOW_Q30:
                    low30 += 1
            lows[mate - 1].append(low)
            lows30[mate - 1].append(low30)
            blocks[mate - 1].append(block_means(q))
        chars = sorted({ch for d in counts.values() for ch in d})
        for b in sorted(counts):
            d = counts[b]
            n = sum(d.values())
            mean = sum((ord(ch) - 33) * k for ch, k in d.items()) / n
            rows.append([mate, b, n, fmt(mean)] + [f"Q{ord(ch) - 33}:{fmt(d.get(ch, 0) / n)}" for ch in chars])
    write_tsv(os.path.join(outdir, "base_quality.tsv"), ["mate", "base", "bases", "mean_q", "quality_fractions"], rows)

    a, b = lows
    if len(a) != len(b) or len(a) < 2:
        raise SystemExit(f"basequal: {r1} and {r2} hold {len(a)} and {len(b)} reads")
    n = len(a)
    corr = statistics.correlation(a, b) if statistics.pstdev(a) > 0 and statistics.pstdev(b) > 0 else 0.0
    a30, b30 = lows30
    corr30 = statistics.correlation(a30, b30) if statistics.pstdev(a30) > 0 and statistics.pstdev(b30) > 0 else 0.0
    block_corr = sum(pearson(x, y) for x, y in zip(*blocks)) / n
    poor1 = sum(1 for x in a if x >= POOR_PAIR) / n
    poor2 = sum(1 for x in b if x >= POOR_PAIR) / n
    both = sum(1 for x, y in zip(a, b) if x >= POOR_PAIR and y >= POOR_PAIR) / n
    write_tsv(
        os.path.join(outdir, "pair_coupling.tsv"),
        ["key", "value"],
        [
            ["pairs", n],
            ["lowq_threshold", LOW_Q],
            ["poor_mate_min_lowq_bases", POOR_PAIR],
            ["lowq_count_correlation", fmt(corr)],
            ["frac_r1_poor", fmt(poor1)],
            ["frac_r2_poor", fmt(poor2)],
            ["frac_both_poor", fmt(both)],
            ["frac_both_poor_if_independent", fmt(poor1 * poor2)],
            ["both_poor_ratio", fmt(both / (poor1 * poor2)) if poor1 * poor2 > 0 else ""],
            ["lowq30_count_correlation", fmt(corr30)],
            ["block_mean_correlation", fmt(block_corr)],
        ],
    )


# ---------------------------------------------------------------- classifiers


def comparisons(cs):
    return {c["feature_set"]: c for c in cs or []}


# Feature sets of profile_matrix.tsv and profile_matching.tsv.
MATRIX_SETS = [("1", "sequence"), ("1", "quality"), ("2", "quality"), ("pair", "pair")]
# Statistics of the profile-matching rule: per mate mean quality of the
# reads (fastq_summary.tsv), which a user can compute from their own run.
MATCH_STATS = [("mean_q", 1), ("mean_q", 2)]


def realism_blocks(path):
    """(read, real-vs-real comparisons, real-vs-generated comparisons) of a
    run_e2e_tests report."""
    with open(path, encoding="utf-8") as f:
        r = json.load(f)
    x = r["platforms"][0]["realism"]
    blocks = [("1", x["real_vs_real"], x.get("real_vs_generated"))]
    for key, read in (("read2", "2"), ("pair", "pair")):
        if x.get(key):
            blocks.append((read, x[key]["real_vs_real"], x[key].get("real_vs_generated")))
    return blocks


def classifiers(outdir, args):
    """classifiers OUTDIR STATS_ROOT PROFILE=RUN,... -- PROFILE@RUN=E2E.json ... -- A:B=E2E.json ..."""
    if args.count("--") != 2:
        raise SystemExit(__doc__)
    i = args.index("--")
    j = args.index("--", i + 1)
    root, prof_arg = args[0], args[1]
    profiles = dict(p.split("=", 1) for p in prof_arg.split(","))
    fitted_runs = set(profiles.values())
    gens = [a.split("=", 1) for a in args[i + 1 : j]]
    crosses = [a.split("=", 1) for a in args[j + 1 :]]
    rows = []
    acc = {}  # (read, features, a, b) -> accuracy; a is "real:RUN" or a profile
    spread = {}  # (read, features) -> kind -> [accuracy]
    runs = []

    def add(kind, a, b, read, c):
        rows.append([kind, a, b, read, c["feature_set"], c["classifier"], c["holdout_reads"], fmt(c["accuracy"], 4), fmt(c["chance_band"], 2)])
        acc[(read, c["feature_set"], a, b)] = c["accuracy"]
        spread.setdefault((read, c["feature_set"]), {}).setdefault(kind, []).append(c["accuracy"])

    seen_null = set()
    for key, path in gens:
        prof, run = key.split("@", 1)
        if run not in runs:
            runs.append(run)
        for read, null, gen in realism_blocks(path):
            if run not in seen_null:
                for c in null:
                    add("real_vs_real", "real:" + run, run, read, c)
            if profiles[prof] == run:
                kind = "generated_vs_own_run"
            elif run in fitted_runs:
                kind = "generated_vs_other_fitted_run"
            else:
                kind = "generated_vs_heldout_run"
            for c in gen or []:
                add(kind, prof, run, read, c)
        seen_null.add(run)
    for pair, path in crosses:
        a, b = pair.split(":", 1)
        with open(path, encoding="utf-8") as f:
            x = json.load(f)["cross_run"]
        for key, read in (("read1", "1"), ("read2", "2"), ("pair", "pair")):
            for c in x.get(key) or []:
                kind = "real_run_vs_real_run"
                add(kind, "real:" + a, b, read, c)
                acc[(read, c["feature_set"], "real:" + b, a)] = c["accuracy"]
    write_tsv(
        os.path.join(outdir, "classifier_matrix.tsv"),
        ["comparison", "a", "b", "read", "features", "classifier", "holdout", "accuracy", "chance_band_4sigma"],
        rows,
    )
    kinds = ["real_vs_real", "generated_vs_own_run", "generated_vs_other_fitted_run", "generated_vs_heldout_run", "real_run_vs_real_run"]
    out = []
    for (read, fs), d in spread.items():
        row = [read, fs]
        for k in kinds:
            v = sorted(d.get(k, []))
            row += [len(v), fmt(v[0], 4) if v else "", fmt(statistics.median(v), 4) if v else "", fmt(v[-1], 4) if v else ""]
        out.append(row)
    header = ["read", "features"]
    for k in kinds:
        header += [f"{k}_n", f"{k}_min", f"{k}_median", f"{k}_max"]
    write_tsv(os.path.join(outdir, "classifier_spread.tsv"), header, out)

    # Profile x run matrix: generated reads of each profile, and each real
    # run, against each run (the diagonal of the real rows is the
    # real-vs-real null).
    sources = list(profiles) + ["real:" + r for r in runs]
    mrows = []
    for read, fs in MATRIX_SETS:
        for s in sources:
            role = "profile fitted to " + profiles[s] if s in profiles else ("fitted run" if s[5:] in fitted_runs else "held-out run")
            mrows.append([read, fs, s, role] + [fmt(acc.get((read, fs, s, r)), 4) for r in runs])
    write_tsv(os.path.join(outdir, "profile_matrix.tsv"), ["read", "features", "source", "role"] + runs, mrows)

    # Profile matching: per run, the profile whose fitted run is closest in
    # the per-mate mean quality, against the best profile and the real runs.
    st = {r: kv_stats(os.path.join(root, r)) for r in runs}
    prows = []
    for run in runs:
        dist = {}
        for p, fr in profiles.items():
            dist[p] = math.sqrt(sum((st[run][k] - st[fr][k]) ** 2 for k in MATCH_STATS))
        # a fitted run is matched among the other profiles
        cands = [p for p in profiles if profiles[p] != run]
        matched = min(cands, key=lambda p: dist[p])
        for read, fs in MATRIX_SETS:
            gen = {p: acc.get((read, fs, p, run)) for p in profiles}
            other = [p for p in cands if gen[p] is not None]
            best = min(other, key=lambda p: gen[p]) if other else ""
            reals = sorted(v for (rd, f, a, b), v in acc.items() if rd == read and f == fs and b == run and a.startswith("real:") and a[5:] != run)
            m = gen.get(matched)
            prows.append([
                run, "fitted" if run in fitted_runs else "held-out", read, fs,
                matched, fmt(dist[matched], 3), fmt(m, 4), best, fmt(gen.get(best), 4),
                len(reals), fmt(reals[0], 4) if reals else "", fmt(statistics.median(reals), 4) if reals else "", fmt(reals[-1], 4) if reals else "",
                ("yes" if m <= reals[-1] else "no") if reals and m is not None else "",
                ("yes" if m <= statistics.median(reals) else "no") if reals and m is not None else "",
            ])
    write_tsv(
        os.path.join(outdir, "profile_matching.tsv"),
        ["run", "role", "read", "features", "matched_profile", "mean_q_distance", "matched_accuracy", "best_other_profile", "best_other_accuracy",
         "real_runs_n", "real_runs_min", "real_runs_median", "real_runs_max", "matched_within_real_max", "matched_within_real_median"],
        prows,
    )


def kv_stats(d):
    """The matching statistics (MATCH_STATS) of a read set."""
    out = {}
    for stat, mate in MATCH_STATS:
        out[(stat, mate)] = float(kv(os.path.join(d, "fastq_summary.tsv"), mate)[stat])
    return out


# ---------------------------------------------------------------- calibration


def kv(path, mate=None):
    out = {}
    for r in read_tsv(path):
        if mate is None:
            out[r["key"]] = r["value"]
        elif r["mate"] == str(mate):
            out[r["key"]] = r["value"]
    return out


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def set_stats(d):
    """Statistics of one read set: {(statistic, mate): value}."""
    s = {}
    for mate in (1, 2):
        fq = kv(os.path.join(d, "fastq_summary.tsv"), mate)
        s[("read_length", mate)] = num(fq.get("max_length"))
        s[("alphabet", mate)] = fq.get("alphabet")
        for k in ("mean_q", "frac_ge_q30", "frac_Q2", "frac_Q11", "frac_Q25", "frac_Q37", "n_rate"):
            s[(k, mate)] = num(fq.get(k, 0))
        # transitions: P(Q11 at c+1 | Q11 at c), quality-change fraction
        same11 = from11 = changes = total = 0
        for r in read_tsv(os.path.join(d, "transitions.tsv")):
            if r["mate"] != str(mate):
                continue
            n = int(r["count"])
            total += n
            if r["from"] != r["to"]:
                changes += n
            if r["from"] == "11":
                from11 += n
                if r["to"] == "11":
                    same11 += n
        s[("p_q11_after_q11", mate)] = same11 / from11 if from11 else None
        s[("qual_change_frac", mate)] = changes / total if total else None
        runs = [(int(r["length"]), int(r["runs"])) for r in read_tsv(os.path.join(d, "run_lengths.tsv")) if r["mate"] == str(mate) and r["phred"] == "11"]
        nrun = sum(n for _, n in runs)
        reads = num(fq.get("reads"))
        s[("q11_runs_per_read", mate)] = nrun / reads if reads else None
        s[("q11_run_mean_length", mate)] = sum(ln * n for ln, n in runs) / nrun if nrun else None
        s[("q11_runs_ge10_per_1000_reads", mate)] = 1000 * sum(n for ln, n in runs if ln >= 10) / reads if reads else None
        s[("q11_runs_ge20_per_1000_reads", mate)] = 1000 * sum(n for ln, n in runs if ln >= 20) / reads if reads else None
        for r in read_tsv(os.path.join(d, "base_quality.tsv")):
            if r["mate"] == str(mate) and r["base"] in "ACGT":
                s[(f"mean_q_{r['base']}", mate)] = num(r["mean_q"])
        for r in read_tsv(os.path.join(d, "errors_quality.tsv")):
            if r["mate"] == str(mate):
                s[(f"mismatch_rate_Q{r['phred']}", mate)] = num(r["raw_rate"])
                s[(f"mismatch_rate_minus_variants_Q{r['phred']}", mate)] = num(r["rate_minus_variants"])
                s[(f"overlap_error_rate_Q{r['phred']}", mate)] = num(r["overlap_rate"])
    sam = kv(os.path.join(d, "sam_summary.tsv"))
    for mate in (1, 2):
        for k in ("raw_mismatch_rate", "ins_events_per_base", "del_events_per_base"):
            s[(k, mate)] = num(sam.get(f"mate{mate}_{k}"))
    pairs = num(sam.get("pairs"))
    s[("proper_usable_pair_frac", 0)] = num(sam.get("proper_usable_pairs")) / pairs if pairs else None
    s[("primary_mapped_frac", 0)] = num(sam.get("primary_mapped")) / num(sam.get("primary")) if num(sam.get("primary")) else None
    for k in ("insert_mean", "insert_sd", "insert_median", "insert_p05", "insert_p95", "insert_below_302", "variant_rate", "indel_error_rate_from_overlap"):
        s[(k, 0)] = num(sam.get(k))
    below = s.pop(("insert_below_302", 0))
    s[("insert_frac_below_302", 0)] = below / num(sam.get("insert_pairs")) if below is not None and num(sam.get("insert_pairs")) else None
    pc = kv(os.path.join(d, "pair_coupling.tsv"))
    for k in ("lowq_count_correlation", "lowq30_count_correlation", "frac_both_poor", "both_poor_ratio", "block_mean_correlation"):
        s[(k, 0)] = num(pc.get(k))
    return s


def cycle_fracs(d):
    """{mate: [ {Qx: fraction} per cycle ]} from cycle_quals.tsv."""
    out = {1: [], 2: []}
    for r in read_tsv(os.path.join(d, "cycle_quals.tsv")):
        bases = int(r["bases"])
        out[int(r["mate"])].append({k: int(v) / bases for k, v in r.items() if k.startswith("Q")})
    return out


def cycle_delta(a, b):
    """Mean over cycles both sets have and over quality values of
    |fraction in a - fraction in b|, per mate."""
    res = {}
    for mate in (1, 2):
        ca, cb = a[mate], b[mate]
        n = min(len(ca), len(cb))
        keys = sorted(set().union(*ca[:n], *cb[:n]))
        tot = 0.0
        for c in range(n):
            for k in keys:
                tot += abs(ca[c].get(k, 0) - cb[c].get(k, 0))
        res[mate] = tot / (n * len(keys)) if n else None
    return res


def calibration(outdir, args):
    """calibration OUTDIR STATS_ROOT PROFILE=RUN,... SET ..."""
    root, prof_arg, sets = args[0], args[1], args[2:]
    profiles = dict(p.split("=", 1) for p in prof_arg.split(","))
    fitted = set(profiles.values())
    stats = {n: set_stats(os.path.join(root, n)) for n in sets}
    gen = {p: "generated-" + p for p in profiles if "generated-" + p in stats}
    heldout = [n for n in sets if n not in fitted and not n.startswith("generated-")]
    keys = []
    for n in sets:
        for k in stats[n]:
            if k not in keys:
                keys.append(k)
    rows = []
    for stat, mate in keys:
        vals = [stats[n].get((stat, mate)) for n in sets]
        hv = [stats[n].get((stat, mate)) for n in heldout]
        hv = [v for v in hv if isinstance(v, float)]
        row = [stat, mate if mate else "pair"] + [fmt(v) if not isinstance(v, str) else v for v in vals]
        row += [fmt(min(hv)), fmt(statistics.median(hv)), fmt(max(hv))] if hv else [""] * 3
        for p, g in gen.items():
            a, b = stats[g].get((stat, mate)), stats[profiles[p]].get((stat, mate))
            row.append(fmt(a - b) if isinstance(a, float) and isinstance(b, float) else "")
        rows.append(row)
    write_tsv(
        os.path.join(outdir, "calibration_stats.tsv"),
        ["statistic", "mate"] + sets + ["heldout_min", "heldout_median", "heldout_max"] + [f"{g}_minus_{profiles[p]}" for p, g in gen.items()],
        rows,
    )
    cyc = {n: cycle_fracs(os.path.join(root, n)) for n in sets}
    drows = []
    for i, a in enumerate(sets):
        for b in sets[i + 1 :]:
            d = cycle_delta(cyc[a], cyc[b])
            drows.append([a, b, fmt(d[1]), fmt(d[2])])
    write_tsv(os.path.join(outdir, "cycle_delta.tsv"), ["a", "b", "mean_abs_delta_r1", "mean_abs_delta_r2"], drows)


def main(argv):
    if len(argv) < 3:
        raise SystemExit(__doc__)
    cmd, outdir = argv[1], argv[2]
    os.makedirs(outdir, exist_ok=True)
    if cmd == "basequal" and len(argv) == 5:
        basequal(outdir, argv[3], argv[4])
    elif cmd == "classifiers" and len(argv) > 6:
        classifiers(outdir, argv[3:])
    elif cmd == "calibration" and len(argv) > 5:
        calibration(outdir, argv[3:])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv)
