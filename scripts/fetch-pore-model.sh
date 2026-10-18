#!/usr/bin/env bash
# Download an ONT k-mer pore model (level table) for nanopore signal
# synthesis (internal/genomics/signal).
#
# The tables come from https://github.com/nanoporetech/kmer_models and are
# NOT redistributed with phantom-science. Review ONT's licence
# (LICENSE.md / COPYRIGHT in that repository) before using or sharing them.
#
# Usage:
#   scripts/fetch-pore-model.sh [MODEL] [DEST_DIR]
#   scripts/fetch-pore-model.sh --list
#
#   MODEL     dna_r10.4.1_e8.2_400bps (default) | dna_r10.4.1_e8.2_260bps | rna004
#   DEST_DIR  default: ${PHANTOM_PORE_MODEL_DIR:-$HOME/.cache/phantom-science/pore-models}
#
# Environment:
#   KMER_MODELS_REF   ref (branch, tag or commit) of nanoporetech/kmer_models
#                     (default: master)
#
# The file is saved as DEST_DIR/<MODEL>/9mer_levels_v1.txt and its path is
# printed on stdout. The tables have two columns (kmer, level) in
# standardised units; load them with the matching preset's scaling:
#   p, _ := signal.Preset(signal.PresetDNAR10)
#   m, err := signal.LoadPoreModel(path, p.ModelOptions()...)
set -euo pipefail

REF="${KMER_MODELS_REF:-master}"
BASE="https://raw.githubusercontent.com/nanoporetech/kmer_models/${REF}"
MODELS=(dna_r10.4.1_e8.2_400bps dna_r10.4.1_e8.2_260bps rna004)

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }

case "${1:-}" in
  -h|--help) usage; exit 0 ;;
  --list) printf '%s\n' "${MODELS[@]}"; exit 0 ;;
esac

MODEL="${1:-dna_r10.4.1_e8.2_400bps}"
DEST="${2:-${PHANTOM_PORE_MODEL_DIR:-$HOME/.cache/phantom-science/pore-models}}"

known=0
for m in "${MODELS[@]}"; do
  if [[ "$m" == "$MODEL" ]]; then known=1; fi
done
if [[ $known -eq 0 ]]; then
  echo "error: unknown model '$MODEL' (known: ${MODELS[*]})" >&2
  exit 2
fi

if command -v curl >/dev/null 2>&1; then
  fetch() { curl -fsSL --retry 3 -o "$2" "$1"; }
elif command -v wget >/dev/null 2>&1; then
  fetch() { wget -q -O "$2" "$1"; }
else
  echo "error: need curl or wget" >&2
  exit 2
fi

FILE="9mer_levels_v1.txt"
OUTDIR="$DEST/$MODEL"
OUT="$OUTDIR/$FILE"
mkdir -p "$OUTDIR"
TMP="$(mktemp "$OUTDIR/.${FILE}.XXXXXX")"
trap 'rm -f "$TMP"' EXIT

URL="$BASE/$MODEL/$FILE"
echo "Downloading $URL" >&2
if ! fetch "$URL" "$TMP"; then
  echo "error: download failed: $URL" >&2
  exit 1
fi

# A 9-mer table has 4^9 = 262144 rows of "<9-mer> <level>".
rows=$(grep -c . "$TMP" || true)
bad=$(awk 'NF < 2 || $1 !~ /^[ACGTU][ACGTU][ACGTU][ACGTU][ACGTU][ACGTU][ACGTU][ACGTU][ACGTU]$/ {n++} END {print n+0}' "$TMP")
if [[ "$rows" -ne 262144 || "$bad" -ne 0 ]]; then
  echo "error: unexpected table ($rows rows, $bad malformed; want 262144 9-mer rows)" >&2
  exit 1
fi
mv "$TMP" "$OUT"
trap - EXIT

if command -v shasum >/dev/null 2>&1; then
  sum=$(shasum -a 256 "$OUT" | awk '{print $1}')
elif command -v sha256sum >/dev/null 2>&1; then
  sum=$(sha256sum "$OUT" | awk '{print $1}')
else
  sum="(no sha256 tool)"
fi

cat >&2 <<NOTE
Saved $OUT
  rows: $rows  sha256: $sum  source ref: $REF

NOTE: this pore model is published by Oxford Nanopore Technologies at
https://github.com/nanoporetech/kmer_models and is not part of
phantom-science. Review ONT's licence (LICENSE.md and COPYRIGHT in that
repository) before using or redistributing it, and do not commit it to
this repository.
NOTE
echo "$OUT"
