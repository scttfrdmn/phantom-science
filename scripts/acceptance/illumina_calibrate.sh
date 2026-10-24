#!/usr/bin/env bash
# Calibration (#27): measure real NovaSeq 6000 reads and fit the Illumina
# quality and error model of internal/genomics/illumina to them.
#
# Usage:
#   scripts/acceptance/illumina_calibrate.sh REFERENCE.fa [READS_1.fastq READS_2.fastq]
#
# REFERENCE is GRCh38 (hg38), FASTA; minimap2 indexes it on the fly (about
# 4 minutes and 10 GB of memory). Without read files the first PAIRS read
# pairs of ENA run ACCESSION are downloaded. Steps:
#   1. quality alphabet, per-cycle quality frequencies, transitions, run
#      lengths, per-read composition and N calls of both mates;
#   2. minimap2 -ax sr alignment; substitution rate per quality and cycle,
#      substitution spectrum, mate-overlap error and variant rates, indel
#      rates and insert sizes of proper pairs;
#   3. the read quality classes (RESULTS/fit.tsv);
#   4. qualities by called-base context (two bases before, the base, one
#      after) per 10-cycle block, and within each read class against the
#      calls expected without a context effect; class posteriors by read GC
#      and longest homopolymer (#35); the joint read classes of the two
#      mates of a pair and the measured mate coupling (#35);
#   5. the model fit, written to RESULTS/fit.tsv and to GO_OUT (gofmt'd),
#      with the measured quality by base and context in
#      RESULTS/context_summary.tsv.
#
# Environment (defaults in brackets):
#   ACCESSION [DRR276626]  PAIRS [200000]  THREADS [8]
#   RESULTS   [scripts/acceptance/calibrate/illumina/results]
#   GO_OUT    [internal/genomics/illumina/calibration.go]
#   WORK      [a new temporary directory, removed at exit unless KEEP_WORK=1]
# Run from the repository root. Needs curl and gzip (download only),
# minimap2, samtools and python3 (standard library only).
set -euo pipefail

die() { echo "illumina_calibrate.sh: $*" >&2; exit 1; }

[[ $# -eq 1 || $# -eq 3 ]] || die "usage: $0 REFERENCE.fa [READS_1.fastq READS_2.fastq]"
ref="$1"
[[ -s "$ref" ]] || die "reference $ref not found"
for tool in minimap2 samtools python3; do
	command -v "$tool" >/dev/null 2>&1 || die "$tool not found on PATH"
done

accession="${ACCESSION:-DRR276626}"
pairs="${PAIRS:-200000}"
threads="${THREADS:-8}"
results="${RESULTS:-scripts/acceptance/calibrate/illumina/results}"
go_out="${GO_OUT:-internal/genomics/illumina/calibration.go}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
stats="$here/calibrate/illumina/illumina_stats.py"
fit="$here/calibrate/illumina/fit_illumina.py"

work="${WORK:-}"
made_work=""
if [[ -z "$work" ]]; then
	work="$(mktemp -d)"
	made_work=1
fi
mkdir -p "$work" "$results"
cleanup() {
	if [[ -n "$made_work" && "${KEEP_WORK:-0}" != 1 ]]; then
		rm -rf "$work"
	fi
}
trap cleanup EXIT

ena_api="https://www.ebi.ac.uk/ena/portal/api/filereport"
fields="run_accession,instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,study_title"

if [[ $# -eq 3 ]]; then
	r1="$2"
	r2="$3"
	[[ -s "$r1" && -s "$r2" ]] || die "read files $r1 / $r2 missing or empty"
	source="local files $(basename "$r1"), $(basename "$r2") (stated to be $accession)"
else
	for tool in curl gzip; do
		command -v "$tool" >/dev/null 2>&1 || die "$tool not found on PATH"
	done
	report="$work/ena.tsv"
	curl -fsS "$ena_api?accession=$accession&result=read_run&fields=$fields&format=tsv" -o "$report" ||
		die "ENA file report for $accession failed"
	urls="$(awk -F'\t' 'NR == 2 { print $8 }' "$report")"
	[[ "$urls" == *";"* ]] || die "$accession has no paired FASTQ files at ENA"
	r1="$work/${accession}_1.fastq"
	r2="$work/${accession}_2.fastq"
	lines=$((4 * pairs))
	for i in 1 2; do
		url="https://$(cut -d';' -f"$i" <<<"$urls")"
		out="$work/${accession}_$i.fastq"
		echo "Downloading the first $pairs reads of $url ..."
		# head closes the pipe early; curl and gzip then fail with SIGPIPE,
		# so the line count below is the check.
		{ curl -fsS "$url" | gzip -dc | head -n "$lines" >"$out"; } || true
		[[ "$(wc -l <"$out" | tr -d ' ')" -eq "$lines" ]] || die "download of $url gave fewer than $pairs reads"
	done
	source="ENA $accession, first $pairs read pairs ($(awk -F'\t' 'NR == 2 { print $2 ", " $3 " " $4 ", nominal insert " $5 }' "$report"))"
fi

echo "Measuring qualities ..."
python3 -I "$stats" fastq "$results" "$r1" "$r2"

echo "Aligning to $ref with minimap2 -ax sr (indexing on the fly) ..."
minimap2 -ax sr --MD -t "$threads" "$ref" "$r1" "$r2" 2>"$work/minimap2.log" |
	samtools view -b -o "$work/aln.bam" -
echo "Measuring errors and insert sizes ..."
samtools view "$work/aln.bam" | python3 -I "$stats" sam "$results"

echo "Fitting the read classes ..."
rm -f "$results/context_classes.tsv" "$results/pair_classes.tsv"
python3 -I "$fit" "$results"
echo "Measuring qualities by sequence context ..."
python3 -I "$stats" context "$results" "$r1" "$r2" "$results/fit.tsv"
echo "Measuring how the mates' read classes go together ..."
python3 -I "$stats" pairs "$results" "$r1" "$r2" "$results/fit.tsv"
echo "Fitting the model ..."
python3 -I "$fit" "$results" "$go_out"
if command -v gofmt >/dev/null 2>&1; then
	gofmt -w "$go_out"
fi

{
	printf 'key\tvalue\n'
	printf 'source\t%s\n' "$source"
	printf 'reads_per_mate\t%s\n' "$(($(wc -l <"$r1") / 4))"
	printf 'reads_md5\t%s\n' "$(python3 -I -c 'import hashlib, sys; print(",".join(hashlib.md5(open(f, "rb").read()).hexdigest() for f in sys.argv[1:]))' "$r1" "$r2")"
	printf 'reference\t%s\n' "$(basename "$ref")"
	if [[ -s "$ref.fai" ]]; then
		printf 'reference_sequences\t%s\n' "$(wc -l <"$ref.fai" | tr -d ' ')"
	fi
	printf 'aligner\tminimap2 %s (-ax sr --MD)\n' "$(minimap2 --version)"
	printf 'samtools\t%s\n' "$(samtools --version | head -n 1)"
	printf 'python\t%s\n' "$(python3 --version 2>&1)"
	printf 'date\t%s\n' "$(date -u +%Y-%m-%d)"
} >"$results/run.tsv"

echo "Wrote $results and $go_out"
column -t -s $'\t' "$results/fastq_summary.tsv"
column -t -s $'\t' "$results/errors_quality.tsv"
