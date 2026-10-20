//go:build !release

package clinical

import "fmt"

// CatalogEntry is one built-in HPO-to-variant mapping. Coordinates are
// GRCh38, 1-based, with VCF-normalised (left-aligned, anchor-padded) REF and
// ALT as given by ClinVar's VCF representation.
type CatalogEntry struct {
	Condition   string
	Terms       []string // HPO term IDs; the entry matches when all are present
	HPOLabels   []string // HPO term names, parallel to Terms
	Gene        string
	HGVS        string // ClinVar preferred name
	ClinVarID   string // ClinVar variation (VCV) accession
	Build       string // genome build of Chrom/Pos/Ref/Alt
	Chrom       string
	Pos         int64
	Ref         string
	Alt         string
	Consequence string // Sequence Ontology term, as reported by ClinVar
	ACMGClass   string // ClinVar germline classification
}

// hpoCatalog maps sets of HPO terms to well-established pathogenic variants.
//
// Source, retrieved 2026-10-19:
//   - HPO terms: hp.obo data-version hp/releases/2026-09-01
//     (https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo),
//     cross-checked with https://ontology.jax.org/api/hp/terms/<ID>.
//   - Variants: ClinVar VCV records via NCBI E-utilities
//     (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id=<id>&rettype=vcv&is_variationid
//     and esummary.fcgi?db=clinvar&id=<id>), GRCh38 positionVCF /
//     referenceAlleleVCF / alternateAlleleVCF, molecular consequence and
//     germline classification; see https://www.ncbi.nlm.nih.gov/clinvar/variation/<id>/.
var hpoCatalog = []CatalogEntry{
	// GLUT1 deficiency syndrome: seizures with low CSF glucose.
	// VCV000016113: Pathogenic (criteria provided, multiple submitters, no
	// conflicts); https://www.ncbi.nlm.nih.gov/clinvar/variation/16113/
	{
		Condition: "GLUT1 deficiency syndrome",
		Terms:     []string{"HP:0001250", "HP:0011972"},
		HPOLabels: []string{"Seizure", "Hypoglycorrhachia"},
		Gene:      "SLC2A1", HGVS: "NM_006516.4(SLC2A1):c.940G>A (p.Gly314Ser)",
		ClinVarID: "VCV000016113", Build: "GRCh38",
		Chrom: "chr1", Pos: 42929242, Ref: "C", Alt: "T",
		Consequence: "missense_variant", ACMGClass: "pathogenic",
	},
	// Phenylketonuria: PAH R408W.
	// VCV000000577: Pathogenic (reviewed by expert panel);
	// https://www.ncbi.nlm.nih.gov/clinvar/variation/577/
	{
		Condition: "Phenylketonuria",
		Terms:     []string{"HP:0004923"},
		HPOLabels: []string{"Hyperphenylalaninemia"},
		Gene:      "PAH", HGVS: "NM_000277.3(PAH):c.1222C>T (p.Arg408Trp)",
		ClinVarID: "VCV000000577", Build: "GRCh38",
		Chrom: "chr12", Pos: 102840493, Ref: "G", Alt: "A",
		Consequence: "missense_variant", ACMGClass: "pathogenic",
	},
	// Cystic fibrosis: CFTR F508del, an in-frame 3-bp deletion.
	// VCV000007105: Pathogenic (practice guideline);
	// https://www.ncbi.nlm.nih.gov/clinvar/variation/7105/
	{
		Condition: "Cystic fibrosis",
		Terms:     []string{"HP:0012236"},
		HPOLabels: []string{"Elevated sweat chloride"},
		Gene:      "CFTR", HGVS: "NM_000492.3(CFTR):c.1521_1523del (p.Phe508del)",
		ClinVarID: "VCV000007105", Build: "GRCh38",
		Chrom: "chr7", Pos: 117559590, Ref: "ATCT", Alt: "A",
		Consequence: "inframe_deletion", ACMGClass: "pathogenic",
	},
}

// Catalog returns a copy of the built-in HPO catalogue.
func Catalog() []CatalogEntry {
	out := make([]CatalogEntry, len(hpoCatalog))
	copy(out, hpoCatalog)
	return out
}

// ResolveHPO resolves HPO terms to a plantable heterozygous variant by
// looking them up in the built-in catalogue (deterministic, no I/O). An entry
// matches when ALL of its Terms are present in the input; the first matching
// entry wins. Terms with no catalogue entry are an error: the scenario must
// then specify planted_variant explicitly.
func ResolveHPO(hpoTerms []string) (*PlantedVariantSpec, *CatalogEntry, error) {
	hpoSet := make(map[string]bool, len(hpoTerms))
	for _, t := range hpoTerms {
		hpoSet[t] = true
	}

	for i := range hpoCatalog {
		e := hpoCatalog[i]
		if allPresent(e.Terms, hpoSet) {
			return &PlantedVariantSpec{
				Chrom:       e.Chrom,
				Pos:         e.Pos,
				Ref:         e.Ref,
				Alt:         e.Alt,
				Zygosity:    "heterozygous",
				Consequence: e.Consequence,
				ACMGClass:   e.ACMGClass,
			}, &e, nil
		}
	}

	return nil, nil, fmt.Errorf("no catalog entry found for HPO terms %v; "+
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
