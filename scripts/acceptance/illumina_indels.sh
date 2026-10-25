#!/usr/bin/env bash
# Indel errors (#30): measure the insertion and deletion error rates of
# real NovaSeq 6000 reads where the two mates of a pair overlap, fit the
# indel model of internal/genomics/illumina to them, and (closed loop)
# measure reads generated with the model the same way.
#
# Usage:
#   scripts/acceptance/illumina_indels.sh REFERENCE.fa
#
# REFERENCE is hg38 primary chromosomes (chr1-22, X, Y), uncompressed FASTA
# with a .fai index (illumina_heldout.sh shows how to make it from UCSC
# hg38.fa.gz). Method (illumina_stats.py indels): in a proper pair's
# overlap both mates read the same molecule, so an indel of the sample (or
# of the library) is in both and an indel in one mate only, after
# left-normalisation, is a sequencing error. Errors are counted per
# homopolymer run, per short tandem repeat and per base, with exact Poisson
# 95 % intervals (indel_*.tsv). Steps:
#   1. the first PAIRS read pairs of each run of RUNS and HELDOUT are
#      streamed from ENA (NovaSeq 6000 binned qualities checked, as in
#      illumina_calibrate.sh), aligned with minimap2 -ax sr --MD and
#      measured into RESULTS/RUN;
#   2. the pool: the first POOL_PAIRS pairs of each run of POOL, measured
#      together into RESULTS/pooled (the runs of novaseq6000-pooled);
#   3. fit_indels.py fits the model's shape to the pool and each profile's
#      rates to its run (PROFILES) and writes GO_OUT and RESULTS/indel_fit.tsv;
#      with RECOVERY, phantom-science is built with that model, PAIRS pairs
#      are generated with each PROFILE=PAIRS of RECOVERY and measured
#      (RESULTS/recovery-PROFILE, seed SEED + 1), and the model is
#      refitted with its insertion and deletion rates multiplied by the
#      run's measured rate over the generated reads' (the measurement loses
#      errors, in long repeats and long indels most);
#   4. the closed loop: phantom-science is built with the new GO_OUT; for
#      each PROFILE=PAIRS of GENERATED, generate (paired, the profile, reads
#      sampled from REFERENCE, seed SEED) writes PAIRS pairs, measured as
#      above into RESULTS/generated-PROFILE; for each of FLOOR the same
#      without indel errors (illumina.errors substitutions), measured into
#      RESULTS/floor-PROFILE: the indels the method finds where there are
#      none (alignment artefacts);
#   5. RESULTS/summary.tsv (the rates of every set), RESULTS/runs.tsv (ENA
#      metadata and md5 of the streamed reads) and RESULTS/run.tsv.
# The runs of the pool have long inserts (DRR394916, DRR276626) and few
# overlapping pairs, so they are measured on more pairs (RUN:PAIRS in
# RUNS) for their own profiles; the pool keeps equal numbers of each.
#
# Environment (defaults in brackets):
#   RUNS       ["DRR394916:4000000 DRR276626:4000000 ERR15600747:1000000"]
#   HELDOUT    ["DRR260182:1000000 SRR30646135:1000000 ERR12737258:1000000"],
#              measured only (no profile is fitted to them)
#   POOL       ["DRR394916 DRR276626 ERR15600747"]  POOL_PAIRS [1000000]
#   PROFILES   ["novaseq6000-pooled=pooled novaseq6000-drr394916=DRR394916
#              novaseq6000-drr276626=DRR276626
#              novaseq6000-err15600747=ERR15600747"]
#   GENERATED  ["novaseq6000-pooled=3000000 novaseq6000-err15600747=2000000";
#              "" to skip the closed loop]
#   FLOOR      ["novaseq6000-err15600747=2000000"; "" to skip]
#   RECOVERY   ["novaseq6000-err15600747=2000000"; "" to fit without it]
#   SEED [30]  THREADS [8]
#   RESULTS    [scripts/acceptance/calibrate/illumina/results/indels]
#   GO_OUT     [internal/genomics/illumina/profile_indels.go]
#   MMI        [WORK/reference.mmi, built with minimap2 -x sr -d]
#   WORK       [a new temporary directory, removed at exit unless KEEP_WORK=1;
#              reads, alignments and the index already in WORK are reused]
# Run from the repository root. Needs go, curl, gzip, minimap2, samtools
# and python3 (standard library only). The reads are about 360 MB per
# 1,000,000 reads; the index 9 GB.
set -euo pipefail

die() { echo "illumina_indels.sh: $*" >&2; exit 1; }

[[ $# -eq 1 ]] || die "usage: $0 REFERENCE.fa"
ref="$1"
[[ -s "$ref" ]] || die "reference $ref not found"
[[ -s "$ref.fai" ]] || samtools faidx "$ref"
for tool in go curl gzip minimap2 samtools python3; do
	command -v "$tool" >/dev/null 2>&1 || die "$tool not found on PATH"
done
[[ -f go.mod ]] || die "run from the repository root"

read -r -a runs <<<"${RUNS:-DRR394916:4000000 DRR276626:4000000 ERR15600747:1000000}"
read -r -a heldout <<<"${HELDOUT-DRR260182:1000000 SRR30646135:1000000 ERR12737258:1000000}"
read -r -a pool <<<"${POOL:-DRR394916 DRR276626 ERR15600747}"
pool_pairs="${POOL_PAIRS:-1000000}"
read -r -a profiles <<<"${PROFILES:-novaseq6000-pooled=pooled novaseq6000-drr394916=DRR394916 novaseq6000-drr276626=DRR276626 novaseq6000-err15600747=ERR15600747}"
read -r -a generated <<<"${GENERATED-novaseq6000-pooled=3000000 novaseq6000-err15600747=2000000}"
read -r -a floor <<<"${FLOOR-novaseq6000-err15600747=2000000}"
read -r -a recovery <<<"${RECOVERY-novaseq6000-err15600747=2000000}"
seed="${SEED:-30}"
threads="${THREADS:-8}"
results="${RESULTS:-scripts/acceptance/calibrate/illumina/results/indels}"
go_out="${GO_OUT:-internal/genomics/illumina/profile_indels.go}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
stats="$here/calibrate/illumina/illumina_stats.py"
fit="$here/calibrate/illumina/fit_indels.py"

work="${WORK:-}"
made_work=""
if [[ -z "$work" ]]; then
	work="$(mktemp -d)"
	made_work=1
fi
mkdir -p "$work/reads" "$work/aln" "$work/gen" "$results"
cleanup() {
	if [[ -n "$made_work" && "${KEEP_WORK:-0}" != 1 ]]; then
		rm -rf "$work"
	fi
}
trap cleanup EXIT

mmi="${MMI:-$work/reference.mmi}"
if [[ ! -s "$mmi" ]]; then
	echo "Indexing $ref for minimap2 -x sr ..."
	minimap2 -x sr -d "$mmi" "$ref" 2>"$work/minimap2-index.log"
fi

md5of() {
	if command -v md5sum >/dev/null 2>&1; then
		md5sum "$1" | cut -d' ' -f1
	else
		md5 -q "$1"
	fi
}

ena_api="https://www.ebi.ac.uk/ena/portal/api/filereport"
fields="run_accession,instrument_model,library_layout,library_strategy,study_accession,center_name,library_construction_protocol,nominal_length,fastq_ftp"
# meta FIELD: the value of FIELD in the ENA report $report.
meta() { awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i } NR == 2 { print $c[k] }' k="$1" "$report"; }

printf 'set\trole\tstudy\tcenter\tlibrary_protocol\tnominal_insert\tpairs\tmd5_r1\tmd5_r2\n' >"$results/runs.tsv"
# fetch RUN PAIRS: the first PAIRS pairs of RUN in $work/reads/RUN.PAIRS_{1,2}.fastq
# (cut from a longer download of RUN already in WORK when there is one).
fetch() {
	local acc=$1 pairs=$2 lines=$((4 * $2)) i out url urls a
	report="$work/reads/$acc.ena.tsv"
	if [[ ! -s "$report" ]]; then
		curl -fsS "$ena_api?accession=$acc&result=read_run&fields=$fields&format=tsv" -o "$report" ||
			die "ENA file report for $acc failed"
	fi
	[[ "$(meta instrument_model)" == "Illumina NovaSeq 6000" ]] || die "$acc is $(meta instrument_model), not Illumina NovaSeq 6000"
	[[ "$(meta library_layout)" == PAIRED ]] || die "$acc is not paired"
	urls="$(meta fastq_ftp)"
	[[ "$urls" == *";"* ]] || die "$acc has no paired FASTQ files at ENA"
	for i in 1 2; do
		out="$work/reads/$acc.${pairs}_$i.fastq"
		if [[ -s "$out" && "$(wc -l <"$out" | tr -d ' ')" -eq "$lines" ]]; then
			continue
		fi
		for longer in "$work/reads/$acc".*_"$i".fastq; do
			if [[ -s "$longer" && "$longer" != "$out" && "$(wc -l <"$longer" | tr -d ' ')" -gt "$lines" ]]; then
				head -n "$lines" "$longer" >"$out"
				break
			fi
		done
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
	done
}
# align NAME R1 R2: $work/aln/NAME.bam (reused when present)
align() {
	if [[ ! -s "$work/aln/$1.bam" ]]; then
		echo "Aligning $1 ..."
		minimap2 -ax sr --MD -t "$threads" "$mmi" "$2" "$3" 2>"$work/aln/$1.minimap2.log" |
			samtools view -b -o "$work/aln/$1.bam.tmp" -
		mv "$work/aln/$1.bam.tmp" "$work/aln/$1.bam"
	fi
}
# measure NAME BAM...: illumina_stats.py indels into RESULTS/NAME
measure() {
	local name=$1
	shift
	echo "Measuring $name ..."
	mkdir -p "$results/$name"
	for b in "$@"; do samtools view "$b"; done | python3 -I "$stats" indels "$results/$name"
}
# record SET ROLE ACC PAIRS
record() {
	report="$work/reads/$3.ena.tsv"
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$(meta study_accession)" "$(meta center_name)" \
		"$(meta library_construction_protocol)" "$(meta nominal_length)" "$4" \
		"$(md5of "$work/reads/$3.${4}_1.fastq")" "$(md5of "$work/reads/$3.${4}_2.fastq")" >>"$results/runs.tsv"
}

# 1. Runs.
for spec in "${runs[@]}" ${heldout[@]+"${heldout[@]}"}; do
	acc="${spec%%:*}"
	pairs="${spec#*:}"
	[[ "$spec" == *:* ]] || pairs="$pool_pairs"
	role=fitted
	[[ " ${runs[*]} " == *" $spec "* ]] || role=heldout
	fetch "$acc" "$pairs"
	align "$acc.$pairs" "$work/reads/$acc.${pairs}_1.fastq" "$work/reads/$acc.${pairs}_2.fastq"
	measure "$acc" "$work/aln/$acc.$pairs.bam"
	record "$acc" "$role" "$acc" "$pairs"
done

# 2. The pool.
bams=()
for acc in "${pool[@]}"; do
	fetch "$acc" "$pool_pairs"
	align "$acc.$pool_pairs" "$work/reads/$acc.${pool_pairs}_1.fastq" "$work/reads/$acc.${pool_pairs}_2.fastq"
	bams+=("$work/aln/$acc.$pool_pairs.bam")
	record pooled "pool" "$acc" "$pool_pairs"
done
measure pooled "${bams[@]}"

# simulate NAME PROFILE PAIRS ERRORS SEED: generate, align and measure into RESULTS/NAME.
simulate() {
	local name=$1 p=$2 n=$3 errs=$4 seed=$5
	# Keyed by the binary (the model and the code), so a refit regenerates.
	local out="$work/gen/$p.$n.$errs.$seed.$(md5of "$work/phantom-science" | cut -c1-8)"
	if [[ ! -s "$out/fastq/gen_S1_L001_R2_001.fastq" ]]; then
		mkdir -p "$out"
		printf '{"sample_name":"gen","lane":1,"read_length":151,"reads_per_batch":%s,"total_batches":1,"paired_end":true}\n' "$n" >"$out/config.json"
		echo "Generating $n pairs with $p (errors $errs) ..."
		"$work/phantom-science" generate --config "$out/config.json" --reference "$ref" --illumina-profile "$p" \
			--illumina-errors "$errs" --seed "$seed" --output "$out/fastq" >"$out/generate.log" 2>&1 ||
			die "generate $p failed; see $out/generate.log"
	fi
	align "$(basename "$out")" "$out/fastq/gen_S1_L001_R1_001.fastq" "$out/fastq/gen_S1_L001_R2_001.fastq"
	measure "$name" "$work/aln/$(basename "$out").bam"
	printf '%s\t%s\t\t\t\t\t%s\t%s\t%s\n' "$name" "generated, errors $errs (seed $seed)" "$n" \
		"$(md5of "$out/fastq/gen_S1_L001_R1_001.fastq")" "$(md5of "$out/fastq/gen_S1_L001_R2_001.fastq")" >>"$results/runs.tsv"
}
build() {
	echo "Building phantom-science ..."
	CGO_ENABLED=0 go build -o "$work/phantom-science" .
}

# 3. Fit (and the recovery refit).
checks=()
for spec in ${heldout[@]+"${heldout[@]}"}; do checks+=("@${spec%%:*}"); done
python3 -I "$fit" "$results" "$go_out" pooled "${profiles[@]}" ${checks[@]+"${checks[@]}"}
gofmt -w "$go_out"
if [[ ${#recovery[@]} -gt 0 ]]; then
	build
	recov=()
	for g in "${recovery[@]}"; do
		simulate "recovery-${g%%=*}" "${g%%=*}" "${g#*=}" quality $((seed + 1))
		recov+=("%recovery-${g%%=*}=${g%%=*}")
	done
	python3 -I "$fit" "$results" "$go_out" pooled "${profiles[@]}" ${checks[@]+"${checks[@]}"} "${recov[@]}"
	gofmt -w "$go_out"
fi

if [[ ${#generated[@]} -gt 0 || ${#floor[@]} -gt 0 ]]; then
	build
fi
# 4. Closed loop (and the method's floor).
for g in ${generated[@]+"${generated[@]}"}; do
	simulate "generated-${g%%=*}" "${g%%=*}" "${g#*=}" quality "$seed"
done
for g in ${floor[@]+"${floor[@]}"}; do
	simulate "floor-${g%%=*}" "${g%%=*}" "${g#*=}" substitutions "$seed"
done

# 5. Summary.
python3 -I - "$results" <<'EOF'
import os, sys
res = sys.argv[1]
keys = ["overlap_pairs", "mate_bases", "ins_rate_per_base", "del_rate_per_base", "indel_rate_per_base",
        "indel_length1_rate_per_base", "indel_length2_rate_per_base", "indel_length3+_rate_per_base", "concordant_rate_per_position"]
with open(os.path.join(res, "summary.tsv"), "w", encoding="ascii") as out:
    out.write("set\toverlap_pairs\tmate_bases\tins_events\tins_rate\tins_low\tins_high\tdel_events\tdel_rate\tdel_low\tdel_high\t"
              "indel_rate\tindel_low\tindel_high\tlength1_rate\tlength2_rate\tlength3plus_rate\tsample_indel_rate\n")
    for name in sorted(os.listdir(res)):
        p = os.path.join(res, name, "indel_summary.tsv")
        if not os.path.exists(p):
            continue
        with open(p, encoding="ascii") as f:
            next(f)
            rows = {r[0]: r for r in (l.rstrip("\n").split("\t") for l in f)}
        v = lambda k, i: rows[k][i]
        out.write("\t".join([name, v("overlap_pairs", 2), v("mate_bases", 2),
                             v("ins_rate_per_base", 1), v("ins_rate_per_base", 2), v("ins_rate_per_base", 3), v("ins_rate_per_base", 4),
                             v("del_rate_per_base", 1), v("del_rate_per_base", 2), v("del_rate_per_base", 3), v("del_rate_per_base", 4),
                             v("indel_rate_per_base", 2), v("indel_rate_per_base", 3), v("indel_rate_per_base", 4),
                             v("indel_length1_rate_per_base", 2), v("indel_length2_rate_per_base", 2), v("indel_length3+_rate_per_base", 2),
                             v("concordant_rate_per_position", 2)]) + "\n")
EOF
{
	printf 'key\tvalue\n'
	printf 'runs\t%s\n' "${runs[*]}"
	printf 'heldout\t%s\n' "${heldout[*]-}"
	printf 'pool\t%s (%s pairs each)\n' "${pool[*]}" "$pool_pairs"
	printf 'profiles\t%s\n' "${profiles[*]}"
	printf 'generated\t%s (seed %s)\n' "${generated[*]-}" "$seed"
	printf 'floor\t%s (seed %s, errors substitutions)\n' "${floor[*]-}" "$seed"
	printf 'recovery\t%s (seed %s, the model before the recovery refit)\n' "${recovery[*]-}" "$((seed + 1))"
	printf 'reference\t%s\n' "$(basename "$ref")"
	printf 'reference_sequences\t%s\n' "$(wc -l <"$ref.fai" | tr -d ' ')"
	printf 'aligner\tminimap2 %s (-ax sr --MD)\n' "$(minimap2 --version)"
	printf 'samtools\t%s\n' "$(samtools --version | head -n 1)"
	printf 'go\t%s\n' "$(go version)"
	printf 'commit\t%s\n' "$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
	printf 'python\t%s\n' "$(python3 --version 2>&1)"
	printf 'date\t%s\n' "$(date -u +%Y-%m-%d)"
} >"$results/run.tsv"

echo "Wrote $results and $go_out"
column -t -s $'\t' "$results/summary.tsv"
