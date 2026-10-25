#!/usr/bin/env bash
# Calibration (#27): measure real NovaSeq 6000 reads and fit the Illumina
# quality and error model of internal/genomics/illumina to them: one run
# profile (#35, illumina.profile) per run.
#
# Usage:
#   scripts/acceptance/illumina_calibrate.sh REFERENCE.fa [READS_1.fastq READS_2.fastq]
#
# REFERENCE is GRCh38 (hg38), FASTA; minimap2 indexes it on the fly (about
# 4 minutes and 10 GB of memory) unless MMI names a prebuilt index
# (minimap2 -x sr -d). Without read files the first PAIRS read pairs of ENA
# run ACCESSION are downloaded. The run's ENA metadata (instrument, centre,
# library protocol, study) describe the profile (RESULTS/profile.tsv); a
# run that is not an Illumina NovaSeq 6000 run, or whose first 10,000
# quality lines use other characters than the binned '#', ',', ':', 'F'
# (SRA-Lite runs have '?' and '$'), is rejected.
#
# A pooled profile (#36): ACCESSION may list several runs (space-separated).
# The first PAIRS read pairs of each are downloaded (read files already in
# WORK with that many reads are reused; local read files cannot be given),
# checked as above and concatenated, so the fit sees equal numbers of pairs
# from every run and its read classes span the runs. Each run's ENA
# metadata goes to RESULTS/profile_runs.tsv; profile.tsv lists the runs.
# Steps:
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
#      mates of a pair and the measured mate coupling (#35); the quality
#      structure along the reads (#36, RESULTS/structure.tsv);
#   5. the model fit, written to RESULTS/fit.tsv and to GO_OUT (gofmt'd,
#      the run's Profile), with the measured quality by base and context in
#      RESULTS/context_summary.tsv.
#
# Environment (defaults in brackets; acc is ACCESSION in lower case, or
# "pooled" when ACCESSION lists several runs):
#   ACCESSION [DRR276626]  PAIRS [200000, per run]  THREADS [8]
#   PROFILE   [novaseq6000-acc], the profile's name
#   RESULTS   [scripts/acceptance/calibrate/illumina/results for DRR276626,
#             scripts/acceptance/calibrate/illumina/results/profiles/pooled
#             for several runs, else
#             scripts/acceptance/calibrate/illumina/results/profiles/ACCESSION]
#   GO_OUT    [internal/genomics/illumina/profile_acc.go]
#   MMI       [none: index REFERENCE on the fly]
#   WORK      [a new temporary directory, removed at exit unless KEEP_WORK=1]
# A new profile must be added to the profiles list in
# internal/genomics/illumina/profiles.go. Run from the repository root.
# Needs curl and gzip (ENA metadata and download), minimap2, samtools and
# python3 (standard library only).
set -euo pipefail

die() { echo "illumina_calibrate.sh: $*" >&2; exit 1; }

[[ $# -eq 1 || $# -eq 3 ]] || die "usage: $0 REFERENCE.fa [READS_1.fastq READS_2.fastq]"
ref="$1"
[[ -s "$ref" ]] || die "reference $ref not found"
for tool in curl gzip minimap2 samtools python3; do
	command -v "$tool" >/dev/null 2>&1 || die "$tool not found on PATH"
done

read -r -a accessions <<<"${ACCESSION:-DRR276626}"
[[ ${#accessions[@]} -ge 1 ]] || die "ACCESSION is empty"
pooled=""
if [[ ${#accessions[@]} -gt 1 ]]; then
	pooled=1
	[[ $# -eq 1 ]] || die "local read files cannot be given for a pooled profile (several ACCESSION runs)"
	accession="$(IFS=,; echo "${accessions[*]}")"
	acc=pooled
else
	accession="${accessions[0]}"
	acc="$(tr '[:upper:]' '[:lower:]' <<<"$accession")"
fi
pairs="${PAIRS:-200000}"
threads="${THREADS:-8}"
profile="${PROFILE:-novaseq6000-$acc}"
default_results=scripts/acceptance/calibrate/illumina/results
if [[ -n "$pooled" ]]; then
	default_results="$default_results/profiles/pooled"
elif [[ "$accession" != DRR276626 ]]; then
	default_results="$default_results/profiles/$accession"
fi
results="${RESULTS:-$default_results}"
go_out="${GO_OUT:-internal/genomics/illumina/profile_$acc.go}"
target="${MMI:-$ref}"
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
fields="run_accession,instrument_model,library_layout,library_strategy,nominal_length,read_count,base_count,fastq_ftp,study_title,study_accession,center_name,library_construction_protocol"
# meta FIELD: the value of FIELD in the ENA report $report.
meta() { awk -F'\t' -v k="$1" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i } NR == 2 { print $c[k] }' "$report"; }
# library_of: the library protocol of the run of $report, as submitted.
library_of() {
	local l
	l="$(meta library_construction_protocol)"
	if [[ -z "$l" ]]; then
		l="not recorded"
		if [[ -n "$(meta nominal_length)" ]]; then
			l="$l (nominal insert $(meta nominal_length))"
		fi
	fi
	printf '%s' "$l"
}
# check_alphabet FASTQ: reject a run without NovaSeq binned qualities.
check_alphabet() {
	local a
	a="$(awk 'NR % 4 == 0 && NR <= 40000' "$1" | fold -w1 | sort -u | tr -d '\n')"
	[[ "$a" =~ ^[#,:F]+$ ]] || die "$1: quality alphabet '$a' is not NovaSeq binned qualities (#,:F); SRA-Lite?"
}
lines=$((4 * pairs))
studies=()
if [[ -n "$pooled" ]]; then
	printf 'accession\tstudy\tinstrument\tcenter\tlibrary\tnominal_insert\tread_length\tpairs\tmd5_r1\tmd5_r2\n' >"$results/profile_runs.tsv"
fi
for a in "${accessions[@]}"; do
	report="$work/ena_$a.tsv"
	curl -fsS "$ena_api?accession=$a&result=read_run&fields=$fields&format=tsv" -o "$report" ||
		die "ENA file report for $a failed"
	[[ "$(meta instrument_model)" == "Illumina NovaSeq 6000" ]] || die "$a is $(meta instrument_model), not Illumina NovaSeq 6000"
	[[ "$(meta library_layout)" == PAIRED ]] || die "$a is not paired"
	studies+=("$(meta study_accession)")
	if [[ $# -eq 3 ]]; then
		continue
	fi
	urls="$(meta fastq_ftp)"
	[[ "$urls" == *";"* ]] || die "$a has no paired FASTQ files at ENA"
	for i in 1 2; do
		out="$work/${a}_$i.fastq"
		if [[ -s "$out" && "$(wc -l <"$out" | tr -d ' ')" -eq "$lines" ]]; then
			echo "Reusing $out"
			continue
		fi
		url="https://$(cut -d';' -f"$i" <<<"$urls")"
		echo "Downloading the first $pairs reads of $url ..."
		# head closes the pipe early; curl and gzip then fail with SIGPIPE,
		# so the line count below is the check.
		{ curl -fsS "$url" | gzip -dc | head -n "$lines" >"$out"; } || true
		[[ "$(wc -l <"$out" | tr -d ' ')" -eq "$lines" ]] || die "download of $url gave fewer than $pairs reads"
	done
	check_alphabet "$work/${a}_1.fastq"
	check_alphabet "$work/${a}_2.fastq"
	if [[ -n "$pooled" ]]; then
		len="$(awk 'NR % 4 == 2 { print length($0) }' "$work/${a}_1.fastq" | sort -n | uniq | tr '\n' ',' | sed 's/,$//')"
		md5s="$(python3 -I -c 'import hashlib, sys; print("\t".join(hashlib.md5(open(f, "rb").read()).hexdigest() for f in sys.argv[1:]))' "$work/${a}_1.fastq" "$work/${a}_2.fastq")"
		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$a" "$(meta study_accession)" "$(meta instrument_model)" "$(meta center_name)" \
			"$(library_of)" "$(meta nominal_length)" "$len" "$pairs" "$md5s" >>"$results/profile_runs.tsv"
	fi
done

if [[ $# -eq 3 ]]; then
	r1="$2"
	r2="$3"
	[[ -s "$r1" && -s "$r2" ]] || die "read files $r1 / $r2 missing or empty"
	check_alphabet "$r1"
	check_alphabet "$r2"
	source="local files $(basename "$r1"), $(basename "$r2") (stated to be $accession)"
elif [[ -n "$pooled" ]]; then
	r1="$work/pooled_1.fastq"
	r2="$work/pooled_2.fastq"
	for i in 1 2; do
		files=()
		for a in "${accessions[@]}"; do files+=("$work/${a}_$i.fastq"); done
		cat "${files[@]}" >"$work/pooled_$i.fastq"
	done
	source="ENA ${accessions[*]}, the first $pairs read pairs of each, concatenated (Illumina NovaSeq 6000, PAIRED WGS; see profile_runs.tsv)"
else
	r1="$work/${accession}_1.fastq"
	r2="$work/${accession}_2.fastq"
	source="ENA $accession, first $pairs read pairs ($(meta instrument_model), $(meta library_layout) $(meta library_strategy), nominal insert $(meta nominal_length))"
fi
{
	printf 'key\tvalue\n'
	printf 'name\t%s\n' "$profile"
	printf 'accession\t%s\n' "$accession"
	printf 'study\t%s\n' "$(IFS=,; echo "${studies[*]}")"
	printf 'instrument\t%s\n' "$(meta instrument_model)"
	if [[ -n "$pooled" ]]; then
		printf 'center\t%s\n' "several (profile_runs.tsv)"
		printf 'library\t%s\n' "several (profile_runs.tsv)"
		printf 'id\t%s\n' "$acc"
		printf 'pairs_per_run\t%s\n' "$pairs"
	else
		printf 'center\t%s\n' "$(meta center_name)"
		printf 'library\t%s\n' "$(library_of)"
	fi
	printf 'pairs\t%s\n' "$(($(wc -l <"$r1") / 4))"
	printf 'results\t%s\n' "$results"
} >"$results/profile.tsv"

echo "Measuring qualities ..."
python3 -I "$stats" fastq "$results" "$r1" "$r2"

echo "Aligning to $target with minimap2 -ax sr ..."
minimap2 -ax sr --MD -t "$threads" "$target" "$r1" "$r2" 2>"$work/minimap2.log" |
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
echo "Measuring the quality structure ..."
python3 -I "$stats" structure "$results" "$r1" "$r2"
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
