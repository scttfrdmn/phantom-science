//go:build !release

package clinical

import (
	"bytes"
	"fmt"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// CheckPlantedRef checks that the REF allele of each planted variant matches
// the reference FASTA at its 1-based position. The FASTA may be plain (with
// or without a samtools .fai next to it; without one the file is scanned) or
// gzip-compressed (scanned). Chromosome names are tried as given and with
// the "chr" prefix added or removed, so a GRCh38 "chr7" variant can be
// checked against an Ensembl-style "7" FASTA. Bases are compared ignoring
// case (soft-masking). The first mismatch, missing chromosome or
// out-of-range position is returned as an error.
func CheckPlantedRef(fastaPath string, variants ...*PlantedVariantSpec) error {
	fp, err := reference.NewFastaParser(fastaPath)
	if err != nil {
		return fmt.Errorf("reference %s: %w", fastaPath, err)
	}
	index, err := fp.LoadIndex()
	if err != nil {
		return fmt.Errorf("reference %s: %w", fastaPath, err)
	}
	for _, v := range variants {
		if v == nil {
			continue
		}
		if v.Ref == "" || v.Pos < 1 {
			return fmt.Errorf("planted variant %s:%d has no REF allele or position to check", v.Chrom, v.Pos)
		}
		chrom, ok := matchChromosome(index, v.Chrom)
		if !ok {
			return fmt.Errorf("planted variant %s:%d %s>%s: chromosome %q is not in reference %s",
				v.Chrom, v.Pos, v.Ref, v.Alt, v.Chrom, fastaPath)
		}
		start := v.Pos - 1
		end := start + int64(len(v.Ref))
		got, err := fp.GetSequence(index, chrom, start, end)
		if err != nil {
			return fmt.Errorf("planted variant %s:%d %s>%s: %w", v.Chrom, v.Pos, v.Ref, v.Alt, err)
		}
		if !bytes.EqualFold(got, []byte(v.Ref)) {
			return fmt.Errorf("planted variant %s:%d %s>%s: REF %q does not match reference %s, which has %q at %s:%d-%d (wrong build or coordinates?)",
				v.Chrom, v.Pos, v.Ref, v.Alt, v.Ref, fastaPath, strings.ToUpper(string(got)), chrom, v.Pos, end)
		}
	}
	return nil
}

// matchChromosome returns the reference's name for chrom, trying it as
// given and with the "chr" prefix toggled.
func matchChromosome(index *reference.FastaIndex, chrom string) (string, bool) {
	candidates := []string{chrom}
	if s, ok := strings.CutPrefix(chrom, "chr"); ok {
		candidates = append(candidates, s)
	} else {
		candidates = append(candidates, "chr"+chrom)
	}
	if chrom == "chrM" || chrom == "MT" {
		candidates = append(candidates, "chrM", "MT", "M")
	}
	for _, c := range candidates {
		if _, ok := index.Chromosomes[c]; ok {
			return c, true
		}
	}
	return "", false
}
