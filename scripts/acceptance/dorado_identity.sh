#!/usr/bin/env bash
# Acceptance check for synthetic nanopore signal (issue #14): basecall
# synthetic BLOW5 with Dorado and align the calls to the known source
# sequences.
#
# Real R10.4.1 DNA reads reach ~98.5% median identity with Dorado SUP, but
# signal synthesised from ONT's k-mer level tables with squigulator's
# step-function algorithm does not: squigulator itself reaches ~91.5%
# (Dorado 2.1.2, dna_r10.4.1_e8.2_400bps_sup@v5.2.0). The pass criterion is
# therefore relative to that algorithm, not to real data:
#   >= 99.9%                -> FAIL, noise model too clean
#   >= baseline - tolerance -> PASS
#   below                   -> FAIL, dwell / stall / skip / noise parameters are off
# The baseline is measured live with --squigulator (squigulator simulates the
# same truth and its signal is basecalled the same way), or set with
# --min-identity (default 90.5%, i.e. the squigulator figure above minus one
# point).
#
# Optionally compares per-base dwell (from Dorado's move table) and pA level
# distributions with a real pod5 slice, e.g. from ONT open data
# (s3://ont-open-data/UHRR_HG002_2026.06).
#
# Usage:
#   scripts/acceptance/dorado_identity.sh --blow5 reads.blow5 --truth truth.fa \
#       [--model sup] [--real-pod5 real.pod5] [--outdir DIR] [--threads N] \
#       [--device DEV] [--rna] [--squigulator [PATH]] [--tolerance PTS] \
#       [--min-identity PCT]
#
#   --blow5      synthetic signal (BLOW5/SLOW5); converted to pod5 with blue-crab s2p
#   --truth      source sequences (FASTA or FASTQ) whose names are the read IDs
#                or at least the sequences reads were drawn from
#   --model      Dorado model or complex (default: sup)
#   --real-pod5  real pod5 file or directory to compare dwell/level histograms
#   --outdir     working directory (default: ./dorado_identity_out)
#   --threads    minimap2/samtools threads (default: 8)
#   --device     Dorado --device (e.g. cpu, metal, cuda:0)
#   --rna        direct RNA: minimap2 -x splice -uf -k14 instead of map-ont
#   --squigulator  also simulate the truth with squigulator (binary on PATH, or
#                PATH) using -x dna-r10-prom (rna004-prom with --rna) and use
#                its median identity minus --tolerance as the pass threshold
#   --tolerance  points below the squigulator baseline still passing (default 1.0)
#   --min-identity  threshold in percent without --squigulator (default 90.5)
#
# Requires: blue-crab (pip install blue-crab; slow5tools itself has no s2p),
# dorado, minimap2, samtools, awk. slow5tools is optional: when present the
# input is first re-encoded to zlib + svb-zd, because the PyPI blue-crab build
# cannot read zstd-compressed BLOW5. The level
# comparison additionally needs python3 with the pod5 and numpy packages
# (skipped otherwise). Exits 0 on PASS, 1 on FAIL, 2 on usage/tool errors.
set -euo pipefail

BLOW5="" TRUTH="" MODEL="sup" REAL="" OUTDIR="./dorado_identity_out" THREADS=8 DEVICE="" RNA=0
SQUIG="" TOLERANCE=1.0 MIN_IDENTITY=90.5

die() { echo "error: $*" >&2; exit 2; }
usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --blow5) BLOW5="$2"; shift 2 ;;
    --truth) TRUTH="$2"; shift 2 ;;
    --model) MODEL="$2"; shift 2 ;;
    --real-pod5) REAL="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --device) DEVICE="$2"; shift 2 ;;
    --rna) RNA=1; shift ;;
    --squigulator)
      if [[ $# -gt 1 && "$2" != --* ]]; then SQUIG="$2"; shift 2; else SQUIG="squigulator"; shift; fi ;;
    --tolerance) TOLERANCE="$2"; shift 2 ;;
    --min-identity) MIN_IDENTITY="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown argument: $1" ;;
  esac
done

[[ -n "$BLOW5" && -n "$TRUTH" ]] || { usage >&2; die "--blow5 and --truth are required"; }
[[ -f "$BLOW5" ]] || die "no such file: $BLOW5"
[[ -f "$TRUTH" ]] || die "no such file: $TRUTH"
[[ -z "$REAL" || -e "$REAL" ]] || die "no such file: $REAL"

missing=()
for tool in blue-crab dorado minimap2 samtools awk ${SQUIG:+"$SQUIG"}; do
  command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
done
if [[ ${#missing[@]} -gt 0 ]]; then
  echo "error: missing tools: ${missing[*]}" >&2
  echo "  blue-crab: pip install blue-crab (https://github.com/Psy-Fer/blue-crab)" >&2
  echo "  dorado:     https://github.com/nanoporetech/dorado" >&2
  echo "  minimap2:   https://github.com/lh3/minimap2   samtools: https://www.htslib.org" >&2
  exit 2
fi

preset=(-x map-ont)
[[ $RNA -eq 1 ]] && preset=(-x splice -uf -k14)

# score BLOW5 DIR: convert to pod5, basecall with moves, align to the truth
# and print the identity summary. Sets MEDIAN and CALLS_BAM.
score() {
  local blow5=$1 dir=$2
  mkdir -p "$dir/pod5"
  local pod5="$dir/pod5/synthetic.pod5" calls="$dir/calls.bam" fastq="$dir/calls.fastq" paf="$dir/calls.paf"
  local s2p_in="$blow5"
  if command -v slow5tools >/dev/null 2>&1; then
    echo "== slow5tools view -c zlib -s svb-zd (re-encode for blue-crab)" >&2
    s2p_in="$dir/reencoded.blow5"
    slow5tools view "$blow5" -c zlib -s svb-zd -o "$s2p_in"
  fi
  echo "== blue-crab s2p" >&2
  rm -f "$pod5"
  blue-crab s2p "$s2p_in" -o "$pod5"

  local dorado_args=(basecaller --emit-moves "$MODEL" "$dir/pod5")
  [[ -n "$DEVICE" ]] && dorado_args+=(--device "$DEVICE")
  echo "== dorado ${dorado_args[*]}" >&2
  dorado "${dorado_args[@]}" > "$calls"
  samtools fastq -@ "$THREADS" "$calls" > "$fastq" 2>/dev/null

  echo "== minimap2 ${preset[*]}" >&2
  minimap2 -c --secondary=no -t "$THREADS" "${preset[@]}" "$TRUTH" "$fastq" > "$paf" 2>/dev/null

  local ncalls naln mean wmean qcov
  ncalls=$(awk 'NR % 4 == 1' "$fastq" | wc -l | tr -d ' ')
  # PAF: col10 = matches, col11 = alignment block length; identity = 10/11
  # (BLAST identity, gaps included). Keep the primary (tp:A:P) alignment of
  # each read.
  awk '{
    tp = ""
    for (i = 13; i <= NF; i++) if ($i ~ /^tp:A:/) tp = substr($i, 6)
    if (tp == "P" && !($1 in seen)) { seen[$1] = 1; print $10 / $11, $10, $11, $2, $4 - $3 }
  }' "$paf" | sort -g > "$dir/identity.txt"
  naln=$(wc -l < "$dir/identity.txt" | tr -d ' ')
  [[ "$naln" -gt 0 ]] || { echo "error: no reads aligned ($ncalls basecalled)" >&2; exit 1; }

  read -r MEDIAN mean wmean qcov <<<"$(awk '
    { id[NR] = $1; s += $1; m += $2; b += $3; q += $5 / $4 }
    END {
      med = (NR % 2) ? id[(NR + 1) / 2] : (id[NR / 2] + id[NR / 2 + 1]) / 2
      printf "%.5f %.5f %.5f %.4f\n", med, s / NR, m / b, q / NR
    }' "$dir/identity.txt")"
  CALLS_BAM=$calls

  echo
  echo "basecalled reads:     $ncalls"
  echo "aligned (primary):    $naln"
  echo "median identity:      $(awk -v x="$MEDIAN" 'BEGIN{printf "%.2f%%", 100*x}')"
  echo "mean identity:        $(awk -v x="$mean" 'BEGIN{printf "%.2f%%", 100*x}')"
  echo "aggregate identity:   $(awk -v x="$wmean" 'BEGIN{printf "%.2f%%", 100*x}')"
  echo "mean query coverage:  $(awk -v x="$qcov" 'BEGIN{printf "%.1f%%", 100*x}')"
}

echo "== synthetic signal"
score "$BLOW5" "$OUTDIR"
median=$MEDIAN
CALLS=$CALLS_BAM
POD5="$OUTDIR/pod5/synthetic.pod5"

threshold=$(awk -v x="$MIN_IDENTITY" 'BEGIN{printf "%.5f", x / 100}')
if [[ -n "$SQUIG" ]]; then
  echo
  echo "== squigulator baseline (same truth)"
  profile=dna-r10-prom
  [[ $RNA -eq 1 ]] && profile=rna004-prom
  # squigulator reads FASTA; take the first word of each header as the name.
  awk 'NR == 1 { fq = /^@/ }
       fq && NR % 4 == 1 { print ">" substr($1, 2); next }
       fq && NR % 4 == 2 { print; next }
       !fq { if (/^>/) print $1; else print }' "$TRUTH" > "$OUTDIR/squigulator_truth.fa"
  "$SQUIG" -x "$profile" --full-contigs --ont-friendly=yes --seed 1 \
    -o "$OUTDIR/squigulator.blow5" "$OUTDIR/squigulator_truth.fa" >&2
  score "$OUTDIR/squigulator.blow5" "$OUTDIR/squigulator"
  baseline=$MEDIAN
  threshold=$(awk -v b="$baseline" -v t="$TOLERANCE" 'BEGIN{printf "%.5f", b - t / 100}')
  echo "squigulator baseline: $(awk -v x="$baseline" 'BEGIN{printf "%.2f%%", 100*x}')"
fi
echo "pass threshold:       $(awk -v x="$threshold" 'BEGIN{printf "%.2f%%", 100*x}')"

# Per-base dwell (samples) from the move table: mv:B:c,stride,m1,m2,...;
# each 1 starts a base, its dwell is stride × (moves until the next 1).
dwell_hist() {
  samtools view "$1" | awk '{
    for (i = 12; i <= NF; i++) if ($i ~ /^mv:B:c,/) {
      n = split(substr($i, 8), mv, ",")
      stride = mv[1]; run = 0; started = 0
      for (j = 2; j <= n; j++) {
        if (mv[j] == 1) { if (started) { d = stride * run; h[d]++; tot++; sum += d; ss += d * d } started = 1; run = 1 }
        else run++
      }
    }
  } END {
    if (tot == 0) { print "no move tables"; exit }
    m = sum / tot
    printf "bases %d  mean dwell %.2f samples  sd %.2f  cv %.3f\n", tot, m, sqrt(ss / tot - m * m), sqrt(ss / tot - m * m) / m
    for (d in h) printf "%d\t%d\n", d, h[d] | "sort -n > \"/dev/stderr\""
  }'
}

echo
echo "== dwell (synthetic, from Dorado moves; histogram in $OUTDIR/dwell_synthetic.tsv)"
dwell_hist "$CALLS" 2> "$OUTDIR/dwell_synthetic.tsv"

if [[ -n "$REAL" ]]; then
  REALCALLS="$OUTDIR/real_calls.bam"
  real_args=(basecaller --emit-moves "$MODEL" "$REAL")
  [[ -n "$DEVICE" ]] && real_args+=(--device "$DEVICE")
  echo "== dorado on real pod5" >&2
  dorado "${real_args[@]}" > "$REALCALLS"
  echo "== dwell (real; histogram in $OUTDIR/dwell_real.tsv)"
  dwell_hist "$REALCALLS" 2> "$OUTDIR/dwell_real.tsv"

  if command -v python3 >/dev/null 2>&1 && python3 -I -c 'import pod5, numpy' 2>/dev/null; then
    PYDIR="$(mktemp -d)"
    trap 'rm -rf "$PYDIR"' EXIT
    cat > "$PYDIR/levels.py" <<'PY'
import sys
import numpy as np
import pod5

def quantiles(path, limit=2000):
    vals = []
    with pod5.DatasetReader(path) as ds:
        for i, rec in enumerate(ds):
            if i >= limit:
                break
            vals.append(rec.signal_pa.astype(np.float32))
    v = np.concatenate(vals)
    q = np.percentile(v, [1, 5, 25, 50, 75, 95, 99])
    return len(vals), v.mean(), v.std(), q

for label, path in (("synthetic", sys.argv[1]), ("real", sys.argv[2])):
    n, m, s, q = quantiles(path)
    print(f"{label:9s} reads {n:5d}  mean {m:7.2f} pA  sd {s:6.2f}  "
          "p1/5/25/50/75/95/99 " + " ".join(f"{x:.1f}" for x in q))
PY
    echo "== pA level distribution"
    python3 -I "$PYDIR/levels.py" "$POD5" "$REAL"
  else
    echo "(level comparison skipped: needs python3 with pod5 and numpy)"
  fi
fi

echo
pass=$(awk -v x="$median" -v t="$threshold" 'BEGIN{print (x >= t && x < 0.999) ? 1 : 0}')
if awk -v x="$median" 'BEGIN{exit !(x >= 0.999)}'; then
  echo "FAIL: identity ~100%: the noise model is too clean"
elif [[ "$pass" -ne 1 ]]; then
  echo "FAIL: identity below the threshold: dwell, stall, skip or noise parameters are off"
else
  echo "PASS: identity at or above the threshold (real R10.4.1 SUP reads reach ~98.5%)"
fi
[[ "$pass" -eq 1 ]]
