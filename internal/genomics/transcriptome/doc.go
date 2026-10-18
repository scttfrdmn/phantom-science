// Package transcriptome generates synthetic Oxford Nanopore cDNA (PCS kit)
// reads from a transcript table, with expression and isoform ground truth
// (#13).
//
// Inputs are a transcript FASTA (GENCODE transcripts.fa or any FASTA), a
// transcript FASTA plus a GTF (for exon structure), or a genome FASTA plus a
// GTF (transcripts are spliced from the genome; minus-strand transcripts are
// reverse-complemented). Each run:
//
//  1. assigns every transcript a TPM: the default ZipfModel (power-law rank
//     abundance with log-normal scatter) or a loaded profile (ParseProfile);
//  2. optionally derives one novel isoform per transcript (an internal exon
//     skipped, or an alternative 3′ end) that takes a NovelRate share of the
//     transcript's molecules;
//  3. per read: draws a molecule species by abundance (so per-species
//     counts are multinomial in TPM), a 5′ truncation (full length with
//     FullLengthFraction, else a geometric 3′ fragment), a poly(A) length,
//     adds the kit's strand-switch and RT primers, picks a random strand,
//     applies the errmodel tier and draws MinKNOW read headers;
//  4. writes FASTQ (optionally gzip), sequencing_summary.txt, truth.tsv,
//     counts.tsv and novel_isoforms.tsv/.gtf.
//
// All randomness derives from the run seed through independent PCG streams
// (expression, isoform design, reads, headers): the same inputs and seed
// give byte-identical output.
package transcriptome
