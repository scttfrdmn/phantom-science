#!/usr/bin/env bash
# Acceptance check (#12, #16): run PycoQC on a sequencing_summary.txt written
# by `phantom-science generate` for an ONT config, and fail unless PycoQC
# parses it and writes a report.
#
# Usage:
#   scripts/acceptance/pycoqc.sh [summary-or-output-dir] [report-dir]
#
# Without arguments a small ONT run is generated first (needs the
# phantom-science binary on PATH or PHANTOM_BIN, and an ONT config in
# PHANTOM_ONT_CONFIG, default configs/sequencers/ont_promethion24.yaml).
# PycoQC is not installed by this script: `pip install pycoQC`.
set -euo pipefail

die() { echo "pycoqc.sh: $*" >&2; exit 1; }

command -v pycoQC >/dev/null 2>&1 || die "pycoQC not found on PATH; install it with: pip install pycoQC"

input="${1:-}"
report_dir="${2:-}"
work=""
cleanup() { [[ -n "$work" ]] && rm -rf "$work"; }
trap cleanup EXIT

if [[ -z "$input" ]]; then
	bin="${PHANTOM_BIN:-phantom-science}"
	command -v "$bin" >/dev/null 2>&1 || die "$bin not found; build it (make build) or set PHANTOM_BIN"
	config="${PHANTOM_ONT_CONFIG:-configs/sequencers/ont_promethion24.yaml}"
	[[ -f "$config" ]] || die "ONT config $config not found; set PHANTOM_ONT_CONFIG"
	work="$(mktemp -d)"
	echo "Generating a 20 Mb ONT run in $work/out ..."
	"$bin" generate --config "$config" --platform ont --total-bases 20M --workers 4 --seed 16 --output "$work/out" >/dev/null
	input="$work/out"
fi

if [[ -d "$input" ]]; then
	summary="$(find "$input" -name 'sequencing_summary*.txt' -print | sort | head -n 1)"
	[[ -n "$summary" ]] || die "no sequencing_summary*.txt under $input"
else
	summary="$input"
fi
[[ -s "$summary" ]] || die "$summary is missing or empty"

if [[ -z "$report_dir" ]]; then
	report_dir="$(mktemp -d)" # kept: only the generated run is removed
fi
mkdir -p "$report_dir"
html="$report_dir/pycoqc.html"
json="$report_dir/pycoqc.json"

rows=$(($(wc -l <"$summary") - 1))
echo "Running PycoQC on $summary ($rows reads) ..."
if ! pycoQC --summary_file "$summary" --html_outfile "$html" --json_outfile "$json" --quiet; then
	die "PycoQC failed on $summary"
fi
[[ -s "$json" ]] || die "PycoQC wrote no JSON report"

# The report must describe every read of the summary.
python3 - "$json" "$rows" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))
want = int(sys.argv[2])
reads = report.get("All Reads", {}).get("basecall", {}).get("reads_number")
if reads is None:
    sys.exit("pycoqc.sh: PycoQC JSON has no 'All Reads.basecall.reads_number'")
if int(reads) != want:
    sys.exit(f"pycoqc.sh: PycoQC counted {reads} reads, the summary has {want}")
print(f"PycoQC parsed all {reads} reads")
PY

echo "PASS: PycoQC report in $html"
