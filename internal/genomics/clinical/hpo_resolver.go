//go:build !release

package clinical

import "fmt"

type hpoEntry struct {
	Gene        string
	Chrom       string
	Pos         int64
	Ref         string
	Alt         string
	Consequence string
	ACMGClass   string
}

// hpoCatalog maps sets of HPO terms to demonstration variants. An entry
// matches when ALL of its Terms are present in the input set. The term IDs
// and variants have not been checked against HPO or ClinVar, and at least
// two IDs do not mean what the entry's comment intends (noted below).
var hpoCatalog = []struct {
	Terms []string
	Entry hpoEntry
}{
	// Intended as GLUT1 deficiency syndrome: seizures (HP:0001250) plus
	// HP:0001253, labelled here as elevated lactate (not verified).
	{
		[]string{"HP:0001250", "HP:0001253"},
		hpoEntry{
			Gene: "SLC2A1", Chrom: "chr1", Pos: 43393429,
			Ref: "C", Alt: "T", Consequence: "missense_variant", ACMGClass: "pathogenic",
		},
	},
	// Intended as PKU (abnormal phenylalanine), but HP:0001360 is
	// Holoprosencephaly in HPO.
	{
		[]string{"HP:0001360"},
		hpoEntry{
			Gene: "PAH", Chrom: "chr12", Pos: 102836819,
			Ref: "G", Alt: "A", Consequence: "missense_variant", ACMGClass: "pathogenic",
		},
	},
	// Intended as cystic fibrosis; HP:0002837 is Recurrent bronchitis in HPO.
	{
		[]string{"HP:0002837"},
		hpoEntry{
			Gene: "CFTR", Chrom: "chr7", Pos: 117559590,
			Ref: "CTT", Alt: "C", Consequence: "frameshift_variant", ACMGClass: "pathogenic",
		},
	},
}

// ResolveHPO resolves HPO terms to a plantable variant specification by
// looking them up in hpoCatalog (deterministic, no I/O). An entry matches
// when ALL of its Terms are present in the input. Terms with no catalog entry
// are an error: the scenario must specify planted_variant explicitly.
func ResolveHPO(hpoTerms []string) (*PlantedVariantSpec, string, error) {
	hpoSet := make(map[string]bool, len(hpoTerms))
	for _, t := range hpoTerms {
		hpoSet[t] = true
	}

	for _, entry := range hpoCatalog {
		if allPresent(entry.Terms, hpoSet) {
			e := entry.Entry
			return &PlantedVariantSpec{
				Chrom:       e.Chrom,
				Pos:         e.Pos,
				Ref:         e.Ref,
				Alt:         e.Alt,
				Zygosity:    "heterozygous",
				Consequence: e.Consequence,
				ACMGClass:   e.ACMGClass,
			}, e.Gene, nil
		}
	}

	return nil, "", fmt.Errorf("no catalog entry found for HPO terms %v; "+
		"specify planted_variant explicitly in scenario YAML", hpoTerms)
}

// allPresent returns true when every term in required is in present.
func allPresent(required []string, present map[string]bool) bool {
	for _, t := range required {
		if !present[t] {
			return false
		}
	}
	return len(required) > 0
}
