#!/usr/bin/env bash
# Acceptance test for the cDNA transcript sampler (#13): simulate ONT cDNA
# reads from a genome + annotation, run the standard ONT long-read
# transcriptomics pipeline (Pychopper -> minimap2 -ax splice -> Bambu), and
# score Bambu's quantification and novel-isoform discovery against truth.
#
# Run by hand, not in CI. Required tools (checked before anything runs):
#   go (to build phantom-science, unless PHANTOM points at a binary)
#   pychopper  >= 2.7   (pip/conda: pychopper)
#   minimap2   >= 2.22
#   samtools   >= 1.10
#   Rscript with the Bioconductor package bambu (>= 3.0)
#
# Usage:
#   scripts/acceptance/transcriptome_bambu.sh GENOME.fa ANNOTATION.gtf [OUTDIR]
#
# A single chromosome keeps it quick, e.g. GRCh38 chr22 and the GENCODE
# primary-assembly GTF filtered to chr22:
#   awk '$1=="chr22"' gencode.v46.primary_assembly.annotation.gtf > chr22.gtf
#
# Environment overrides (defaults in brackets):
#   READS [200000]  SEED [1]  NOVEL_RATE [0.05]  TIER [sup]  KIT [PCS114]
#   THREADS [8]  NDR [Bambu default]  MIN_NOVEL_READS [5]
#   PHANTOM [build ./phantom-science into OUTDIR/bin]
#
# Metrics written to OUTDIR/acceptance.tsv and printed:
#   - Pearson correlation of log1p(true reads) vs log1p(Bambu counts) and
#     Spearman correlation, over annotated transcripts;
#   - novel-isoform precision (Bambu novel transcripts whose intron chain
#     equals an injected exon-skip isoform) and recall (injected exon-skip
#     isoforms with >= MIN_NOVEL_READS reads that Bambu reports).
# Alternative-3'-end isoforms keep their parent's intron chain, so a
# junction-based caller like Bambu cannot report them as new; they are
# counted but not scored.

set -euo pipefail

usage() {
	sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
	exit 2
}

[[ $# -ge 2 && $# -le 3 ]] || usage
GENOME=$1
GTF=$2
OUT=${3:-transcriptome_bambu_out}
READS=${READS:-200000}
SEED=${SEED:-1}
NOVEL_RATE=${NOVEL_RATE:-0.05}
TIER=${TIER:-sup}
KIT=${KIT:-PCS114}
THREADS=${THREADS:-8}
NDR=${NDR:-}
MIN_NOVEL_READS=${MIN_NOVEL_READS:-5}

die() {
	echo "transcriptome_bambu: $*" >&2
	exit 1
}

[[ -r $GENOME ]] || die "cannot read genome FASTA $GENOME"
[[ -r $GTF ]] || die "cannot read annotation GTF $GTF"

missing=()
for tool in pychopper minimap2 samtools Rscript; do
	command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
done
if [[ -z ${PHANTOM:-} ]]; then
	command -v go >/dev/null 2>&1 || missing+=("go (or set PHANTOM)")
fi
if [[ ${#missing[@]} -gt 0 ]]; then
	die "missing required tools: ${missing[*]}"
fi
Rscript -e 'suppressPackageStartupMessages(library(bambu))' >/dev/null 2>&1 ||
	die "R package bambu is not installed (BiocManager::install(\"bambu\"))"

case $KIT in
PCS109) PYCHOPPER_KIT=PCS109 ;;
PCS110) PYCHOPPER_KIT=PCS110 ;;
PCS111 | PCS114) PYCHOPPER_KIT=PCS111 ;;
*) die "unknown KIT $KIT (want PCS109, PCS110, PCS111 or PCS114)" ;;
esac

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
REPO=$(cd "$(dirname "$0")/../.." && pwd)

if [[ -z ${PHANTOM:-} ]]; then
	echo "== building phantom-science"
	mkdir -p "$OUT/bin"
	(cd "$REPO" && CGO_ENABLED=0 go build -o "$OUT/bin/phantom-science" .)
	PHANTOM=$OUT/bin/phantom-science
fi
[[ -x $PHANTOM ]] || die "PHANTOM=$PHANTOM is not executable"

echo "== simulating $READS cDNA reads (seed $SEED, novel rate $NOVEL_RATE, tier $TIER, kit $KIT)"
"$PHANTOM" transcriptome --genome "$GENOME" --gtf "$GTF" --reads "$READS" \
	--out "$OUT/sim" --seed "$SEED" --novel-rate "$NOVEL_RATE" --tier "$TIER" --kit "$KIT"

echo "== pychopper (-k $PYCHOPPER_KIT)"
pychopper -k "$PYCHOPPER_KIT" -m edlib -t "$THREADS" \
	-r "$OUT/pychopper_report.pdf" -S "$OUT/pychopper_stats.tsv" \
	-u "$OUT/unclassified.fq" -w "$OUT/rescued.fq" \
	"$OUT/sim/reads.fastq" "$OUT/full_length.fq"
cat "$OUT/full_length.fq" "$OUT/rescued.fq" >"$OUT/pychopped.fq"

echo "== minimap2 -ax splice"
minimap2 -ax splice -uf -k14 -t "$THREADS" "$GENOME" "$OUT/pychopped.fq" |
	samtools sort -@ "$THREADS" -o "$OUT/aln.bam" -
samtools index "$OUT/aln.bam"
[[ -e $GENOME.fai ]] || samtools faidx "$GENOME" ||
	die "cannot index $GENOME (Bambu needs $GENOME.fai; make the directory writable or index it)"

echo "== bambu"
mkdir -p "$OUT/bambu"
cat >"$OUT/score_bambu.R" <<'RSCRIPT'
suppressPackageStartupMessages(library(bambu))
a <- commandArgs(trailingOnly = TRUE)
bam <- a[1]; gtf <- a[2]; genome <- a[3]; outdir <- a[4]; simdir <- a[5]
ncore <- as.integer(a[6]); ndr <- a[7]; min_reads <- as.integer(a[8]); metrics <- a[9]

anno <- prepareAnnotations(gtf)
se <- if (nzchar(ndr)) {
  bambu(reads = bam, annotations = anno, genome = genome, ncore = ncore, NDR = as.numeric(ndr))
} else {
  bambu(reads = bam, annotations = anno, genome = genome, ncore = ncore)
}
writeBambuOutput(se, path = outdir)

tx <- rownames(se)
est <- as.numeric(SummarizedExperiment::assays(se)$counts[, 1])
names(est) <- tx

rd <- function(f) read.delim(file.path(simdir, f), quote = "", colClasses = "character", na.strings = NULL)
truth <- rd("counts.tsv")
truth$reads <- as.numeric(truth$reads)
ann <- truth[truth$novel_isoform == "", ]
e <- est[ann$transcript_id]; e[is.na(e)] <- 0
pearson <- cor(log1p(ann$reads), log1p(e))
spearman <- cor(ann$reads, e, method = "spearman")

# Intron chains: "chrom:strand:start-end,..." from Bambu's exon ranges.
rr <- SummarizedExperiment::rowRanges(se)
chain <- function(gr) {
  gr <- sort(gr)
  n <- length(gr)
  if (n < 2) return(NA_character_)
  paste0(as.character(GenomicRanges::seqnames(gr))[1], ":", as.character(GenomicRanges::strand(gr))[1], ":",
         paste0(GenomicRanges::end(gr)[-n] + 1, "-", GenomicRanges::start(gr)[-1] - 1, collapse = ","))
}
mc <- S4Vectors::mcols(rr)
is_novel <- if ("novelTranscript" %in% names(mc)) as.logical(mc$novelTranscript) else grepl("^Bambu", names(rr))
bambu_keys <- vapply(as.list(rr[is_novel]), chain, character(1))
bambu_keys <- unique(bambu_keys[!is.na(bambu_keys)])

nov <- rd("novel_isoforms.tsv")
nov$reads <- as.numeric(nov$reads)
skip <- nov[nov$novel_isoform == "skip_exon" & nov$intron_chain != "", ]
skip$key <- paste0(skip$chrom, ":", skip$strand, ":", skip$intron_chain)
detectable <- skip[skip$reads >= min_reads, ]
recall <- if (nrow(detectable) > 0) mean(detectable$key %in% bambu_keys) else NA
precision <- if (length(bambu_keys) > 0) mean(bambu_keys %in% skip$key) else NA

m <- data.frame(
  metric = c("annotated_transcripts", "pearson_log1p_counts", "spearman_counts",
             "bambu_novel_multiexon", "injected_skip_isoforms", "injected_skip_isoforms_detectable",
             "injected_alt3_isoforms", "novel_precision", "novel_recall"),
  value = c(nrow(ann), pearson, spearman, length(bambu_keys), nrow(skip), nrow(detectable),
            sum(nov$novel_isoform == "alt_3prime"), precision, recall))
write.table(m, metrics, sep = "\t", quote = FALSE, row.names = FALSE)
print(m, row.names = FALSE)
RSCRIPT
Rscript --vanilla "$OUT/score_bambu.R" "$OUT/aln.bam" "$GTF" "$GENOME" "$OUT/bambu" "$OUT/sim" \
	"$THREADS" "$NDR" "$MIN_NOVEL_READS" "$OUT/acceptance.tsv"

echo "== done: metrics in $OUT/acceptance.tsv"
