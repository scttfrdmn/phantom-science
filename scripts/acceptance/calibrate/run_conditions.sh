#!/usr/bin/env bash
# Run signal-calibration conditions (issue #28) from conditions.tsv, one
# after another (Dorado uses the GPU), appending rows to a results TSV.
#
# Usage:
#   run_conditions.sh [--inputs] [--conditions FILE] [--results FILE] [NAME...]
#
#   --inputs      first (re)build $INPUTS/hg38.fa and uniform.fa with
#                 make_inputs.py and check their sha256 against conditions.tsv
#   --conditions  condition table (default: conditions.tsv next to this script)
#   --results     results TSV (default: results/dna.tsv next to this script)
#   NAME...       run only these conditions (default: all; a condition whose
#                 name is already in the results TSV is skipped)
#
# Environment: PS (/tmp/ps-accept), INPUTS ($PS/w2e/inputs), TABLES
# ($PS/w2e/tables, derived pore tables), PHANTOM, and everything eval.sh
# reads. $PS, $INPUTS and $TABLES in the args column are expanded.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
export PS="${PS:-/tmp/ps-accept}"
INPUTS="${INPUTS:-$PS/w2e/inputs}"
TABLES="${TABLES:-$PS/w2e/tables}"
# Shared eval.sh arguments of the RNA004 conditions ($RNAARGS).
RNAARGS="--rna --model $PS/models/rna004_130bps_sup@v5.3.0 --table $PS/models/rna004/9mer_levels_v1.txt --table-scale 17 --table-shift 80 --mm2-arg -x --mm2-arg map-ont --mm2-arg -k5 --mm2-arg -w1 --mm2-arg -m10 --mm2-arg -s15 --mm2-arg -n2"
PHANTOM="${PHANTOM:-$PS/w2e/bin/phantom-science}"
CONDS="$HERE/conditions.tsv" RESULTS="$HERE/results/dna.tsv" BUILD=0
NAMES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --inputs) BUILD=1; shift ;;
    --conditions) CONDS="$2"; shift 2 ;;
    --results) RESULTS="$2"; shift 2 ;;
    -h|--help) awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit 0 ;;
    -*) echo "error: unknown argument $1" >&2; exit 2 ;;
    *) NAMES+=("$1"); shift ;;
  esac
done

if [[ $BUILD -eq 1 ]]; then
  mkdir -p "$INPUTS"
  python3 -I "$HERE/make_inputs.py" hg38 --fasta "$PS/ref/hg38.fa" \
    --regions "$HERE/regions.tsv" --out "$INPUTS/hg38.fa"
  "$PHANTOM" generate --platform ont --config "$REPO/configs/sequencers/ont_promethion24.yaml" \
    --total-bases 11M --seed 1 --stdout > "$INPUTS/generate.fastq" 2> "$INPUTS/generate.log"
  python3 -I "$HERE/make_inputs.py" uniform --like "$INPUTS/hg38.fa" \
    --fastq "$INPUTS/generate.fastq" --out "$INPUTS/uniform.fa"
  rm -f "$INPUTS/generate.fastq"
  for f in hg38.fa uniform.fa; do
    want=$(awk -v f="$f" '$2 == f { s = 1 } s && $2 == "sha256" { print $3; exit }' "$CONDS")
    got=$(shasum -a 256 "$INPUTS/$f" | cut -d' ' -f1)
    [[ -z "$want" || "$want" == "$got" ]] || { echo "error: $f sha256 $got, conditions.tsv pins $want" >&2; exit 1; }
  done
fi

while IFS=$'\t' read -r name sim input args; do
  [[ -z "$name" || "$name" == \#* ]] && continue
  if [[ ${#NAMES[@]} -gt 0 ]]; then
    hit=0
    for n in "${NAMES[@]}"; do [[ "$n" == "$name" ]] && hit=1; done
    [[ $hit -eq 1 ]] || continue
  fi
  if [[ -f "$RESULTS" ]] && awk -F'\t' -v n="$name" '$1 == n { f = 1 } END { exit !f }' "$RESULTS"; then
    echo "== $name: already in $(basename "$RESULTS"), skipped" >&2
    continue
  fi
  args="${args//\$PS/$PS}"
  args="${args//\$INPUTS/$INPUTS}"
  args="${args//\$TABLES/$TABLES}"
  args="${args//\$RNAARGS/$RNAARGS}"
  in="$input"
  [[ "$input" == "-" ]] && in="real" || in="$INPUTS/$input"
  # shellcheck disable=SC2086 # args is a list of eval.sh options
  "$HERE/eval.sh" --condition "$name" --sim "$sim" --input "$in" --results "$RESULTS" $args
done < "$CONDS"
