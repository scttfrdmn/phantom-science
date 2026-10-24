#!/usr/bin/env bash
# Held-out validation (#35): measure each Illumina run profile, fitted to
# one NovaSeq 6000 run, against its own run, the other profiles' runs and
# NovaSeq 6000 runs no profile was fitted to, and against the distance
# between real runs.
#
# Usage:
#   scripts/acceptance/illumina_heldout.sh [REFERENCE.fa]
#
# REFERENCE is hg38 primary chromosomes (chr1-22, X, Y), uncompressed FASTA
# with a .fai index. Without it UCSC hg38.fa.gz is downloaded into WORK and
# the primary chromosomes are extracted with samtools faidx. Steps:
#   1. the first PAIRS read pairs of each run in RUNS are streamed from ENA
#      (the ENA metadata of each run is recorded); a run whose first 10,000
#      quality lines use other characters than the NovaSeq binned alphabet
#      '#', ',', ':', 'F' is rejected (SRA-Lite runs have only '?' and '$');
#   2. cmd/run_e2e_tests -platforms novaseq -reference REFERENCE
#      -illumina-profile P on each run for each profile P, paired
#      (-real-fastq2): real vs real (alternate halves) and real vs the
#      profile's generated reads for read 1, read 2 and pairs; the
#      calibration (null, GC and quality controls, pair set included) for
#      the first profile on its run; each profile's generated reads for its
#      own run are kept (generated-P);
#   3. run_e2e_tests -real-only -real-fastq-b for every pair of runs: real
#      run A vs real run B, the realistic target;
#   4. illumina_stats.py (fastq, and sam on minimap2 -ax sr alignments to
#      REFERENCE) and heldout_compare.py basequal on every run and on every
#      profile's reads;
#   5. heldout_compare.py writes the tables to RESULTS: classifier_matrix,
#      classifier_spread, profile_matrix, profile_matching,
#      calibration_stats, cycle_delta, runs (ENA metadata, alphabet, md5 of
#      the streamed reads) and run.tsv; and for each read set its small
#      tables under RESULTS/sets/NAME.
#
# Environment (defaults in brackets):
#   RUNS     ["DRR276626 DRR394916 ERR15600747 DRR260182 SRR30646135
#            ERR12737258"]; the runs of PROFILES and the held-out runs
#   PROFILES ["novaseq6000-drr276626=DRR276626 novaseq6000-drr394916=DRR394916
#            novaseq6000-err15600747=ERR15600747"]: each profile and the run
#            it was fitted to; runs without a profile are held out
#   PAIRS    [200000]  THREADS [8]  SEED [1]
#   RESULTS [scripts/acceptance/calibrate/illumina/results/heldout]
#   MMI     [WORK/reference.mmi, built with minimap2 -x sr -d]
#   WORK    [a new temporary directory, removed at exit unless KEEP_WORK=1;
#           reads and the index already in WORK are reused]
# Run from the repository root. Needs go, curl, gzip, minimap2, samtools,
# md5 or md5sum, and python3 (standard library only). The downloads are
# about 70 MB per read file (3 GB more for hg38, 9 GB for the index).
set -euo pipefail

die() { echo "illumina_heldout.sh: $*" >&2; exit 1; }

[[ $# -le 1 ]] || die "usage: $0 [REFERENCE.fa]"
for tool in go curl gzip minimap2 samtools python3; do
	command -v "$tool" >/dev/null 2>&1 || die "$tool not found on PATH"
done
[[ -f go.mod ]] || die "run from the repository root"

read -r -a runs <<<"${RUNS:-DRR276626 DRR394916 ERR15600747 DRR260182 SRR30646135 ERR12737258}"
[[ ${#runs[@]} -ge 2 ]] || die "RUNS needs at least two runs"
read -r -a profiles <<<"${PROFILES:-novaseq6000-drr276626=DRR276626 novaseq6000-drr394916=DRR394916 novaseq6000-err15600747=ERR15600747}"
[[ ${#profiles[@]} -ge 1 ]] || die "PROFILES needs at least one profile"
prof_names=()
for pr in "${profiles[@]}"; do
	[[ "$pr" == *=* ]] || die "PROFILES entry $pr is not PROFILE=RUN"
	[[ " ${runs[*]} " == *" ${pr#*=} "* ]] || die "the run of profile $pr is not in RUNS"
	prof_names+=("${pr%%=*}")
done
prof_list="$(IFS=,; echo "${profiles[*]}")"
pairs="${PAIRS:-200000}"
threads="${THREADS:-8}"
seed="${SEED:-1}"
results="${RESULTS:-scripts/acceptance/calibrate/illumina/results/heldout}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
stats="$here/calibrate/illumina/illumina_stats.py"
compare="$here/calibrate/illumina/heldout_compare.py"
first_profile="${prof_names[0]}"

work="${WORK:-}"
made_work=""
if [[ -z "$work" ]]; then
	work="$(mktemp -d)"
	made_work=1
fi
mkdir -p "$work/reads" "$work/e2e" "$work/stats" "$results"
cleanup() {
	if [[ -n "$made_work" && "${KEEP_WORK:-0}" != 1 ]]; then
		rm -rf "$work"
	fi
}
trap cleanup EXIT

md5of() {
	if command -v md5sum >/dev/null 2>&1; then
		md5sum "$1" | cut -d' ' -f1
	else
		md5 -q "$1"
	fi
}

# Reference.
if [[ $# -eq 1 ]]; then
	ref="$1"
	[[ -s "$ref" ]] || die "reference $ref not found"
	[[ -s "$ref.fai" ]] || samtools faidx "$ref"
else
	ref="$work/hg38.primary.fa"
	if [[ ! -s "$ref.fai" ]]; then
		echo "Downloading UCSC hg38 ..."
		curl -fsS -o "$work/hg38.fa.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
		gzip -dc "$work/hg38.fa.gz" >"$work/hg38.fa"
		samtools faidx "$work/hg38.fa"
		chroms=()
		for c in $(seq 1 22) X Y; do chroms+=("chr$c"); done
		samtools faidx "$work/hg38.fa" "${chroms[@]}" >"$ref"
		samtools faidx "$ref"
		rm -f "$work/hg38.fa" "$work/hg38.fa.fai" "$work/hg38.fa.gz"
	fi
fi
mmi="${MMI:-$work/reference.mmi}"
if [[ ! -s "$mmi" ]]; then
	echo "Indexing $ref for minimap2 -x sr ..."
	minimap2 -x sr -d "$mmi" "$ref" 2>"$work/minimap2-index.log"
fi

# 1. Reads.
ena_api="https://www.ebi.ac.uk/ena/portal/api/filereport"
fields="run_accession,instrument_model,library_layout,library_strategy,library_source,study_accession,center_name,library_construction_protocol,nominal_length,read_count,base_count,first_public,fastq_ftp"
lines=$((4 * pairs))
# meta FIELD: the value of FIELD in the ENA report $report.
meta() { awk -F'\t' -v k="$1" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i } NR == 2 { print $c[k] }' "$report"; }
printf 'run\tstudy\tcenter\tinstrument\tlayout\tstrategy\tsource\tlibrary_protocol\tnominal_insert\tread_length\talphabet_r1\talphabet_r2\tpairs\tmd5_r1\tmd5_r2\n' >"$results/runs.tsv"
for acc in "${runs[@]}"; do
	report="$work/reads/$acc.ena.tsv"
	curl -fsS "$ena_api?accession=$acc&result=read_run&fields=$fields&format=tsv" -o "$report" ||
		die "ENA file report for $acc failed"
	[[ "$(meta instrument_model)" == "Illumina NovaSeq 6000" ]] || die "$acc is $(meta instrument_model), not Illumina NovaSeq 6000"
	[[ "$(meta library_layout)" == PAIRED ]] || die "$acc is not paired"
	urls="$(meta fastq_ftp)"
	[[ "$urls" == *";"* ]] || die "$acc has no paired FASTQ files at ENA"
	alpha=()
	for i in 1 2; do
		out="$work/reads/${acc}_$i.fastq"
		if [[ ! -s "$out" || "$(wc -l <"$out" | tr -d ' ')" -ne "$lines" ]]; then
			url="https://$(cut -d';' -f"$i" <<<"$urls")"
			echo "Downloading the first $pairs reads of $url ..."
			# head closes the pipe early; curl and gzip then fail with
			# SIGPIPE, so the line count below is the check.
			{ curl -fsS "$url" | gzip -dc | head -n "$lines" >"$out"; } || true
			[[ "$(wc -l <"$out" | tr -d ' ')" -eq "$lines" ]] || die "download of $url gave fewer than $pairs reads"
		fi
		a="$(awk 'NR % 4 == 0 && NR <= 40000' "$out" | fold -w1 | sort -u | tr -d '\n')"
		[[ "$a" =~ ^[#,:F]+$ ]] || die "$acc read $i: quality alphabet '$a' is not NovaSeq binned qualities (#,:F); SRA-Lite?"
		alpha+=("$a")
	done
	len="$(awk 'NR % 4 == 2 { print length($0) }' "$work/reads/${acc}_1.fastq" | sort -n | uniq | tr '\n' ',' | sed 's/,$//')"
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$acc" "$(meta study_accession)" "$(meta center_name)" \
		"$(meta instrument_model)" "$(meta library_layout)" "$(meta library_strategy)" "$(meta library_source)" \
		"$(meta library_construction_protocol)" "$(meta nominal_length)" "$len" "${alpha[0]}" "${alpha[1]}" "$pairs" \
		"$(md5of "$work/reads/${acc}_1.fastq")" "$(md5of "$work/reads/${acc}_2.fastq")" >>"$results/runs.tsv"
done

# 2-3. Classifiers.
echo "Building run_e2e_tests ..."
CGO_ENABLED=0 go build -o "$work/run_e2e_tests" ./cmd/run_e2e_tests
harness() {
	local name=$1
	shift
	echo "run_e2e_tests: $name ..."
	"$work/run_e2e_tests" -platforms novaseq -reads "$pairs" -seed "$seed" -benchmark=false \
		-reference "$ref" -output "$work/e2e/$name" "$@" >"$work/e2e/$name.txt" 2>"$work/e2e/$name.log" ||
		die "run_e2e_tests $name failed; see $work/e2e/$name.txt"
}
run_args=()
for pr in "${profiles[@]}"; do
	p="${pr%%=*}"
	own="${pr#*=}"
	for acc in "${runs[@]}"; do
		auth=false
		gen=()
		if [[ "$acc" == "$own" ]]; then
			gen=(-write-generated "$work/generated-$p")
			if [[ "$p" == "$first_profile" ]]; then
				auth=true
			fi
		fi
		harness "${p}_$acc" -illumina-profile "$p" -authenticity="$auth" ${gen[@]+"${gen[@]}"} \
			-real-fastq "$work/reads/${acc}_1.fastq" -real-fastq2 "$work/reads/${acc}_2.fastq"
		run_args+=("$p@$acc=$work/e2e/${p}_$acc/e2e_results.json")
	done
done
cross_args=()
for ((i = 0; i < ${#runs[@]}; i++)); do
	for ((j = i + 1; j < ${#runs[@]}; j++)); do
		a="${runs[i]}"
		b="${runs[j]}"
		harness "${a}_vs_$b" -authenticity=false -real-only \
			-real-fastq "$work/reads/${a}_1.fastq" -real-fastq2 "$work/reads/${a}_2.fastq" \
			-real-fastq-b "$work/reads/${b}_1.fastq" -real-fastq2-b "$work/reads/${b}_2.fastq"
		cross_args+=("$a:$b=$work/e2e/${a}_vs_$b/e2e_results.json")
	done
done
for pr in "${profiles[@]}"; do
	p="${pr%%=*}"
	cp "$work/e2e/${p}_${pr#*=}.txt" "$results/run_e2e_tests_${p}_${pr#*=}.txt"
done

# 4. Calibration statistics of every run and of every profile's reads.
models=()
for p in "${prof_names[@]}"; do
	mv "$work/generated-$p/illumina_novaseq_6000_1.fastq" "$work/reads/generated-${p}_1.fastq"
	mv "$work/generated-$p/illumina_novaseq_6000_2.fastq" "$work/reads/generated-${p}_2.fastq"
	models+=("generated-$p")
done
for name in "${models[@]}" "${runs[@]}"; do
	echo "Measuring $name ..."
	out="$work/stats/$name"
	mkdir -p "$out"
	python3 -I "$stats" fastq "$out" "$work/reads/${name}_1.fastq" "$work/reads/${name}_2.fastq"
	python3 -I "$compare" basequal "$out" "$work/reads/${name}_1.fastq" "$work/reads/${name}_2.fastq"
	minimap2 -ax sr --MD -t "$threads" "$mmi" "$work/reads/${name}_1.fastq" "$work/reads/${name}_2.fastq" 2>"$out/minimap2.log" |
		python3 -I "$stats" sam "$out"
	mkdir -p "$results/sets/$name"
	for f in fastq_summary sam_summary errors_quality cycle_quals run_lengths insert_sizes base_quality pair_coupling; do
		cp "$out/$f.tsv" "$results/sets/$name/"
	done
done

# 5. Tables.
python3 -I "$compare" classifiers "$results" "$work/stats" "$prof_list" -- "${run_args[@]}" -- "${cross_args[@]}"
python3 -I "$compare" calibration "$results" "$work/stats" "$prof_list" "${models[@]}" "${runs[@]}"
{
	printf 'key\tvalue\n'
	printf 'runs\t%s\n' "${runs[*]}"
	printf 'profiles\t%s\n' "${profiles[*]}"
	printf 'pairs_per_run\t%s\n' "$pairs"
	printf 'seed\t%s\n' "$seed"
	printf 'model_reads\trun_e2e_tests -platforms novaseq -reference -illumina-profile P -real-fastq RUN (paired, binned model, the profile'"'"'s insert sizes), %s pairs; generated-P against the run of P\n' "$pairs"
	printf 'reference\t%s\n' "$(basename "$ref")"
	printf 'reference_sequences\t%s\n' "$(wc -l <"$ref.fai" | tr -d ' ')"
	printf 'aligner\tminimap2 %s (-ax sr --MD)\n' "$(minimap2 --version)"
	printf 'samtools\t%s\n' "$(samtools --version | head -n 1)"
	printf 'go\t%s\n' "$(go version)"
	printf 'commit\t%s\n' "$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
	printf 'python\t%s\n' "$(python3 --version 2>&1)"
	printf 'date\t%s\n' "$(date -u +%Y-%m-%d)"
} >"$results/run.tsv"

echo "Wrote $results"
column -t -s $'\t' "$results/classifier_spread.tsv"
column -t -s $'\t' "$results/profile_matching.tsv"
