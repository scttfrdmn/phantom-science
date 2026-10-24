"""Held-out validation of the Illumina model (#35).

Run with `python3 -I` (no third-party packages). Subcommands:

  basequal OUTDIR R1.fastq R2.fastq
      Per mate and called base: bases, mean quality and the fraction of
      each quality value (base_quality.tsv); how the two mates' qualities
      go together: the correlation of the mates' low-quality (Phred < 20)
      base counts and the pairs with both mates poor against the number
      expected if the mates were independent (pair_coupling.tsv).

  classifiers OUTDIR NAME=E2E_RESULTS.json ... -- A:B=E2E_RESULTS.json ...
      Reads run_e2e_tests reports: one per run (real vs real and real vs
      generated) and one per pair of runs (-real-fastq-b, real run A vs
      real run B). Writes classifier_matrix.tsv (every comparison) and
      classifier_spread.tsv (per read and feature set: the real-vs-real
      nulls, generated vs each held-out run, generated vs the fitted run,
      and real run vs real run, as min / median / max).

  calibration OUTDIR MODEL FITTED STATS_ROOT NAME ...
      Compares calibration statistics of each read set (illumina_stats.py
      fastq and sam output, and basequal output, in STATS_ROOT/NAME) with
      the generated reads (MODEL) and the fitted run (FITTED). Writes
      calibration_stats.tsv (one row per statistic and mate, one column per
      read set, plus the spread of the held-out runs) and cycle_delta.tsv
      (mean |difference| of the per-cycle quality fractions between every
      two read sets).

All output is tab-separated with a header line.
"""

import json
import math
import os
import statistics
import sys

LOW_Q = 20
POOR_PAIR = 5  # a mate is poor with at least this many bases below LOW_Q


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
    for mate, path in enumerate((r1, r2), start=1):
        counts = {}  # base -> {qual char: n}
        for s, q in fastq_records(path):
            low = 0
            for b, ch in zip(s, q):
                d = counts.setdefault(b, {})
                d[ch] = d.get(ch, 0) + 1
                if ord(ch) - 33 < LOW_Q:
                    low += 1
            lows[mate - 1].append(low)
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
        ],
    )


# ---------------------------------------------------------------- classifiers


def comparisons(cs):
    return {c["feature_set"]: c for c in cs or []}


def classifiers(outdir, args):
    if "--" not in args:
        raise SystemExit(__doc__)
    sep = args.index("--")
    runs = [a.split("=", 1) for a in args[:sep]]
    crosses = [a.split("=", 1) for a in args[sep + 1 :]]
    rows = []
    # (read, set) -> kind -> [accuracy]
    spread = {}
    fitted = runs[0][0]

    def add(kind, a, b, read, c):
        rows.append([kind, a, b, read, c["feature_set"], c["classifier"], c["holdout_reads"], fmt(c["accuracy"], 4), fmt(c["chance_band"], 2)])
        spread.setdefault((read, c["feature_set"]), {}).setdefault(kind, []).append(c["accuracy"])

    for name, path in runs:
        with open(path, encoding="utf-8") as f:
            r = json.load(f)
        x = r["platforms"][0]["realism"]
        blocks = [("1", x["real_vs_real"], x["real_vs_generated"])]
        for key, read in (("read2", "2"), ("pair", "pair")):
            if x.get(key):
                blocks.append((read, x[key]["real_vs_real"], x[key].get("real_vs_generated")))
        for read, null, gen in blocks:
            for c in null:
                add("real_vs_real", name, name, read, c)
            for c in gen or []:
                add("generated_vs_fitted_run" if name == fitted else "generated_vs_heldout_run", name, "generated", read, c)
    for pair, path in crosses:
        a, b = pair.split(":", 1)
        with open(path, encoding="utf-8") as f:
            x = json.load(f)["cross_run"]
        for key, read in (("read1", "1"), ("read2", "2"), ("pair", "pair")):
            for c in x.get(key) or []:
                add("real_run_vs_real_run", a, b, read, c)
    write_tsv(
        os.path.join(outdir, "classifier_matrix.tsv"),
        ["comparison", "a", "b", "read", "features", "classifier", "holdout", "accuracy", "chance_band_4sigma"],
        rows,
    )
    kinds = ["real_vs_real", "generated_vs_fitted_run", "generated_vs_heldout_run", "real_run_vs_real_run"]
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
    for k in ("lowq_count_correlation", "frac_both_poor", "both_poor_ratio"):
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


def calibration(outdir, model, fitted, root, names):
    sets = [model] + names
    stats = {n: set_stats(os.path.join(root, n)) for n in sets}
    heldout = [n for n in names if n != fitted]
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
        if hv:
            row += [fmt(min(hv)), fmt(statistics.median(hv)), fmt(max(hv))]
            m, f = stats[model].get((stat, mate)), stats[fitted].get((stat, mate))
            row += [fmt(m - statistics.median(hv)) if isinstance(m, float) else "", fmt(f - statistics.median(hv)) if isinstance(f, float) else ""]
        else:
            row += [""] * 5
        rows.append(row)
    write_tsv(
        os.path.join(outdir, "calibration_stats.tsv"),
        ["statistic", "mate"] + sets + ["heldout_min", "heldout_median", "heldout_max", "model_minus_heldout_median", "fitted_minus_heldout_median"],
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
    elif cmd == "classifiers" and len(argv) > 4:
        classifiers(outdir, argv[3:])
    elif cmd == "calibration" and len(argv) > 6:
        calibration(outdir, argv[3], argv[4], argv[5], argv[6:])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv)
