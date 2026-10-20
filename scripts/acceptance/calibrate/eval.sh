#!/usr/bin/env bash
# Score one signal-calibration condition (issue #28) and append a row to a
# results TSV.
#
# A condition is a truth read set plus a simulator: phantom-science signal
# (with --set profile overrides), squigulator, or a real pod5 (reference
# rows). The signal is converted to pod5, basecalled with Dorado
# (--emit-moves), aligned to the truth with minimap2, and scored by
# measure.py: median identity with a bootstrap 95% CI over reads, Dorado
# move-table dwell statistics and per-base signal statistics. Intermediates
# are deleted once the row is written (unless --keep).
#
# Usage:
#   eval.sh --condition NAME --sim phantom|squigulator|real --input FASTA \
#       --results FILE [--set NAME=VALUE]... [--pore-model FILE] \
#       [--profile NAME] [--squig-profile NAME] [--squig-arg ARG]... \
#       [--pod5 FILE] [--ref FASTA|MMI] [--model DIR] [--min-mapq N] \
#       [--mm2-arg ARG]... [--max-reads N] [--keep]
#
#   --input     truth FASTA (synthetic) or, for --sim real, a label
#   --ref       alignment target (default: --input); for real reads GRCh38
#   --pod5      real pod5 file (--sim real)
#   --set       phantom-science signal --profile-set NAME=VALUE
#   --squig-arg extra squigulator argument (repeatable)
#   --model     Dorado model directory (default: $MODELS/dna_r10.4.1_e8.2_400bps_sup@v5.2.0)
#   --min-mapq  minimum MAPQ of scored alignments (default 0; 20 for real reads)
#   --max-reads use only the first N truth records (real: the first N pod5 reads)
#   --table     k-mer level table for measure.py's level-table fit, scaled to
#               pA by --table-scale/--table-shift (default: the R10.4.1 9-mer
#               table, 23.506/93.692); --table none skips the fit
#   --rna       direct RNA: signal runs 3'->5' of the calls (measure.py --rna)
#
# Environment (defaults): WORK (/tmp/ps-accept/w2e/work), PS (/tmp/ps-accept),
# PHANTOM ($PS/w2e/bin/phantom-science), DORADO, SQUIGULATOR, BLUECRAB,
# PYTHON (the pod5/numpy venv python), TOOLS_LOCK ($PS/TOOLS.lock),
# DEVICE (metal), THREADS (8), SEED (1).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
PS="${PS:-/tmp/ps-accept}"
WORK="${WORK:-$PS/w2e/work}"
PHANTOM="${PHANTOM:-$PS/w2e/bin/phantom-science}"
DORADO="${DORADO:-$PS/dorado-2.1.2-osx-arm64/bin/dorado}"
SQUIGULATOR="${SQUIGULATOR:-$PS/bin/squigulator}"
BLUECRAB="${BLUECRAB:-$PS/bc-venv/bin/blue-crab}"
PYTHON="${PYTHON:-$PS/bc-venv/bin/python}"
TOOLS_LOCK="${TOOLS_LOCK:-$PS/TOOLS.lock}"
DEVICE="${DEVICE:-metal}"
THREADS="${THREADS:-8}"
SEED="${SEED:-1}"

COND="" SIM="" INPUT="" RESULTS="" POREMODEL="" PROFILE="" SQPROFILE="dna-r10-prom"
POD5="" REF="" MODEL="$PS/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0" MINMAPQ=0
MM2=() MAXREADS=0 KEEP=0 RNA=0
TABLE="$PS/models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt" TSCALE=23.506 TSHIFT=93.692
SETS=() SQARGS=()

die() { echo "error: $*" >&2; exit 2; }
usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --condition) COND="$2"; shift 2 ;;
    --sim) SIM="$2"; shift 2 ;;
    --input) INPUT="$2"; shift 2 ;;
    --results) RESULTS="$2"; shift 2 ;;
    --set) SETS+=("$2"); shift 2 ;;
    --pore-model) POREMODEL="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --squig-profile) SQPROFILE="$2"; shift 2 ;;
    --squig-arg) SQARGS+=("$2"); shift 2 ;;
    --pod5) POD5="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --model) MODEL="$2"; shift 2 ;;
    --min-mapq) MINMAPQ="$2"; shift 2 ;;
    --mm2-arg) MM2+=("$2"); shift 2 ;;
    --max-reads) MAXREADS="$2"; shift 2 ;;
    --table) TABLE="$2"; shift 2 ;;
    --table-scale) TSCALE="$2"; shift 2 ;;
    --table-shift) TSHIFT="$2"; shift 2 ;;
    --rna) RNA=1; shift ;;
    --keep) KEEP=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown argument: $1" ;;
  esac
done
[[ -n "$COND" && -n "$SIM" && -n "$INPUT" && -n "$RESULTS" ]] || { usage >&2; die "--condition, --sim, --input and --results are required"; }
[[ "$COND" =~ ^[A-Za-z0-9_.+=-]+$ ]] || die "condition name must be [A-Za-z0-9_.+=-]"
for t in "$DORADO" "$BLUECRAB" "$PYTHON" minimap2 samtools; do
  command -v "$t" >/dev/null 2>&1 || die "missing tool: $t"
done
REF="${REF:-$INPUT}"

dir="$WORK/$COND"
rm -rf "$dir"
mkdir -p "$dir/pod5"
log() { echo "[$(date +%T)] $COND: $*" >&2; }

truth="$INPUT"
if [[ "$SIM" != real && "$MAXREADS" -gt 0 ]]; then
  truth="$dir/truth_in.fa"
  awk -v n="$MAXREADS" '/^>/ { if (++c > n) exit } { print }' "$INPUT" > "$truth"
  REF="$truth"
fi
if [[ "$SIM" == real ]]; then
  [[ -f "$POD5" ]] || die "--pod5 is required for --sim real"
  input_sha=$(shasum -a 256 "$POD5" | cut -d' ' -f1)
else
  [[ -f "$truth" ]] || die "no such file: $truth"
  input_sha=$(shasum -a 256 "$truth" | cut -d' ' -f1)
fi

settings=""
case "$SIM" in
  phantom)
    [[ -n "$POREMODEL" ]] || die "--pore-model is required for --sim phantom"
    args=(signal --input "$truth" --out "$dir/sig" --pore-model "$POREMODEL" --seed "$SEED" --reads-per-file 1000000)
    [[ -n "$PROFILE" ]] && args+=(--profile "$PROFILE")
    for s in ${SETS[@]+"${SETS[@]}"}; do args+=(--profile-set "$s"); done
    settings="${PROFILE:-default}${SETS[*]:+ ${SETS[*]}}"
    log "phantom-science ${args[*]}"
    "$PHANTOM" "${args[@]}" > "$dir/signal.log" 2>&1
    blow5=$(ls "$dir"/sig/*.blow5)
    [[ $(wc -l <<<"$blow5") -eq 1 ]] || die "expected one BLOW5 file in $dir/sig"
    ;;
  squigulator)
    settings="squigulator -x $SQPROFILE${SQARGS[*]:+ ${SQARGS[*]}}"
    blow5="$dir/squigulator.blow5"
    log "$settings"
    "$SQUIGULATOR" -x "$SQPROFILE" --full-contigs --ont-friendly=yes --seed "$SEED" \
      ${SQARGS[@]+"${SQARGS[@]}"} -o "$blow5" "$truth" > "$dir/squigulator.log" 2>&1
    ;;
  real)
    [[ -f "$POD5" ]] || die "--pod5 is required for --sim real"
    settings="real $(basename "$POD5")"
    ;;
  *) die "unknown --sim $SIM" ;;
esac

if [[ "$SIM" == real && "$MAXREADS" -gt 0 ]]; then
  # The first N reads in file order.
  "$PYTHON" -I -c 'import sys, pod5
with pod5.Reader(sys.argv[1]) as r:
    for i, rec in zip(range(int(sys.argv[2])), r.reads()):
        print(rec.read_id)' "$POD5" "$MAXREADS" > "$dir/ids.txt"
  "$(dirname "$PYTHON")/pod5" filter "$POD5" --ids "$dir/ids.txt" --output "$dir/pod5/input.pod5" > "$dir/filter.log" 2>&1
  settings="$settings first $MAXREADS reads"
elif [[ "$SIM" == real ]]; then
  ln -s "$POD5" "$dir/pod5/input.pod5"
else
  log "blue-crab s2p"
  "$BLUECRAB" s2p "$blow5" -o "$dir/pod5/input.pod5" > "$dir/s2p.log" 2>&1
  rm -rf "$dir/sig" "$blow5"
fi

log "dorado basecaller $(basename "$MODEL")"
t0=$(date +%s)
"$DORADO" basecaller --emit-moves --device "$DEVICE" "$MODEL" "$dir/pod5" > "$dir/calls.bam" 2> "$dir/dorado.log"
basecall_s=$(( $(date +%s) - t0 ))
samtools view "$dir/calls.bam" > "$dir/calls.sam"
samtools fastq "$dir/calls.bam" > "$dir/calls.fq" 2>/dev/null
rm -f "$dir/calls.bam"
[[ ${#MM2[@]} -gt 0 ]] || MM2=(-x map-ont)
log "minimap2 ${MM2[*]}"
minimap2 -c --secondary=no -t "$THREADS" "${MM2[@]}" "$REF" "$dir/calls.fq" > "$dir/calls.paf" 2> "$dir/minimap2.log"

log "measure"
margs=(--sam "$dir/calls.sam" --paf "$dir/calls.paf" --pod5 "$dir/pod5/input.pod5"
  --min-mapq "$MINMAPQ" --seed "$SEED")
[[ "$TABLE" != none ]] && margs+=(--table "$TABLE" --table-scale "$TSCALE" --table-shift "$TSHIFT")
[[ $RNA -eq 1 ]] && margs+=(--rna)
"$PYTHON" -I "$HERE/measure.py" "${margs[@]}" > "$dir/measure.txt" 2> "$dir/measure.log"
cat "$dir/measure.txt" >&2

cols=(condition input input_sha256 sim settings n_reads n_aligned median_identity ci_lo ci_hi
  mean_identity aligned_frac median_tident bases mean_dwell dwell_cv dwell_short noise_sd
  noise_mad diff_ac1 seg_mean_sd step_abs resid_ac1 carry1 carry2 signal_mean signal_sd
  table_r table_resid_sd table_offset basecall_s dorado_model device git_sha tools_lock_sha256 date)
git_sha=$(git -C "$REPO" rev-parse --short HEAD)
[[ -z "$(git -C "$REPO" status --porcelain -- internal/genomics/signal cmd/signal.go 2>/dev/null)" ]] || git_sha="$git_sha-dirty"
lock_sha=$(shasum -a 256 "$TOOLS_LOCK" | cut -d' ' -f1)
{
  echo "condition=$COND"
  echo "input=$(basename "$INPUT")"
  echo "input_sha256=$input_sha"
  echo "sim=$SIM"
  echo "settings=$settings"
  echo "basecall_s=$basecall_s"
  echo "dorado_model=$(basename "$MODEL")"
  echo "device=$DEVICE"
  echo "git_sha=$git_sha"
  echo "tools_lock_sha256=$lock_sha"
  echo "date=$(date -u +%FT%TZ)"
  cat "$dir/measure.txt"
} > "$dir/row.kv"
mkdir -p "$(dirname "$RESULTS")"
[[ -s "$RESULTS" ]] || (IFS=$'\t'; echo "${cols[*]}") > "$RESULTS"
awk -F= -v cols="${cols[*]}" '
  { k = $1; sub(/^[^=]*=/, ""); v[k] = $0 }
  END { n = split(cols, c, " "); for (i = 1; i <= n; i++) printf "%s%s", (c[i] in v ? v[c[i]] : "NA"), (i < n ? "\t" : "\n") }
' "$dir/row.kv" >> "$RESULTS"
log "row appended to $RESULTS"

if [[ $KEEP -eq 0 ]]; then
  rm -rf "$dir/pod5" "$dir/calls.sam" "$dir/calls.fq"
fi
