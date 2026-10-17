//go:build !release

package clinical

import (
    "context"
    "fmt"

    "github.com/scttfrdmn/phantom-science/internal/genomics/diseases"
)

type hpoEntry struct {
    Gene        string
    Chrom       string
    Pos         int64
    Ref         string
    Alt         string
    Consequence string
    ACMGClass   string
}

// hpoCatalog maps sets of HPO terms to well-known pathogenic variants.
// An entry matches when ALL of its Terms are present in the input set.
var hpoCatalog = []struct {
    Terms []string
    Entry hpoEntry
}{
    // GLUT1 deficiency syndrome: seizures (HP:0001250) + elevated lactate (HP:0001253)
    {
        []string{"HP:0001250", "HP:0001253"},
        hpoEntry{
            Gene: "SLC2A1", Chrom: "chr1", Pos: 43393429,
            Ref: "C", Alt: "T", Consequence: "missense_variant", ACMGClass: "pathogenic",
        },
    },
    // PKU: abnormal phenylalanine (HP:0001360)
    {
        []string{"HP:0001360"},
        hpoEntry{
            Gene: "PAH", Chrom: "chr12", Pos: 102836819,
            Ref: "G", Alt: "A", Consequence: "missense_variant", ACMGClass: "pathogenic",
        },
    },
    // Cystic fibrosis: recurrent pulmonary infections (HP:0002837)
    {
        []string{"HP:0002837"},
        hpoEntry{
            Gene: "CFTR", Chrom: "chr7", Pos: 117559590,
            Ref: "CTT", Alt: "C", Consequence: "frameshift_variant", ACMGClass: "pathogenic",
        },
    },
}

// ResolveHPO resolves HPO terms to a plantable variant specification.
// Resolution order:
//  1. hpoCatalog — deterministic, no I/O; matches when ALL catalog Terms ⊆ input.
//  2. WellCharacterizedDiseases — if a disease's GeneSymbols match HPO term keywords,
//     call FetchPathogenicVariants and convert the first result.
//  3. error with guidance to specify planted_variant explicitly.
func ResolveHPO(ctx context.Context, hpoTerms []string, cv *diseases.ClinVarIntegrator) (*PlantedVariantSpec, string, error) {
    hpoSet := make(map[string]bool, len(hpoTerms))
    for _, t := range hpoTerms {
        hpoSet[t] = true
    }

    // 1. Catalog lookup.
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

    // 2. WellCharacterizedDiseases fallback (requires ClinVar integrator).
    if cv != nil {
        for _, disease := range diseases.WellCharacterizedDiseases {
            if !hpoTermsMatchDisease(hpoTerms, disease.GeneSymbols) {
                continue
            }
            variants, err := cv.FetchPathogenicVariants(ctx, disease.GeneSymbols)
            if err != nil || len(variants) == 0 {
                continue
            }
            v := variants[0]
            if v.Chromosome == "" || v.Position == 0 || v.ReferenceAllele == "" || v.AlternateAllele == "" {
                continue
            }
            return &PlantedVariantSpec{
                Chrom:    v.Chromosome,
                Pos:      v.Position,
                Ref:      v.ReferenceAllele,
                Alt:      v.AlternateAllele,
                Zygosity: "heterozygous",
            }, v.GeneSymbol, nil
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

// hpoTermsMatchDisease returns true if any HPO term string matches a gene symbol
// as a rough keyword heuristic (e.g. "CFTR" in term text).
func hpoTermsMatchDisease(hpoTerms []string, geneSymbols []string) bool {
    for _, term := range hpoTerms {
        for _, gene := range geneSymbols {
            if len(gene) >= 3 && containsInsensitive(term, gene) {
                return true
            }
        }
    }
    return false
}

func containsInsensitive(s, substr string) bool {
    if len(s) < len(substr) {
        return false
    }
    sl := toLower(s)
    subl := toLower(substr)
    for i := 0; i <= len(sl)-len(subl); i++ {
        if sl[i:i+len(subl)] == subl {
            return true
        }
    }
    return false
}

func toLower(s string) string {
    b := make([]byte, len(s))
    for i := range s {
        c := s[i]
        if c >= 'A' && c <= 'Z' {
            c += 'a' - 'A'
        }
        b[i] = c
    }
    return string(b)
}
