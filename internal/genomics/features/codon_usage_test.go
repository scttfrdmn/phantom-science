//go:build !release
// +build !release

package features

import (
	"testing"
)

// Tests for codon usage loading functions

func TestNewCodonUsageModel_NotNil(t *testing.T) {
	model, err := NewCodonUsageModel()
	if err != nil {
		t.Fatalf("NewCodonUsageModel() error: %v", err)
	}
	if model == nil {
		t.Error("NewCodonUsageModel() returned nil")
	}
}

func TestNewCodonUsageModel_HasFrequencies(t *testing.T) {
	model, err := NewCodonUsageModel()
	if err != nil {
		t.Fatalf("NewCodonUsageModel() error: %v", err)
	}
	// NewCodonUsageModel populates expectedFrequencies from human genome data
	if len(model.expectedFrequencies) == 0 {
		t.Error("NewCodonUsageModel() has no expected frequencies")
	}
}

func TestLoadDefaultCodonUsage_Human(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "human"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(human) produced no frequencies")
	}
}

func TestLoadDefaultCodonUsage_HomoSapiens(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "homo_sapiens"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(homo_sapiens) produced no frequencies")
	}
}

func TestLoadDefaultCodonUsage_Mouse(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "mouse"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(mouse) produced no frequencies")
	}
}

func TestLoadDefaultCodonUsage_MusMusculus(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "mus_musculus"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(mus_musculus) produced no frequencies")
	}
}

func TestLoadDefaultCodonUsage_Yeast(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "yeast"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(yeast) produced no frequencies")
	}
}

func TestLoadDefaultCodonUsage_SaccharomycesCerevisiae(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "saccharomyces_cerevisiae"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(saccharomyces_cerevisiae) produced no frequencies")
	}
}

func TestLoadDefaultCodonUsage_EColi(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "ecoli"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(ecoli) produced no frequencies")
	}
}

func TestLoadDefaultCodonUsage_EscherichiaColi(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "escherichia_coli"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(escherichia_coli) produced no frequencies")
	}
}

func TestLoadDefaultCodonUsage_Unknown(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.organism = "unknown_organism"
	model.codonFrequencies = make(map[string]float64)
	model.loadDefaultCodonUsage()
	// Universal codon usage should still populate frequencies
	if len(model.codonFrequencies) == 0 {
		t.Error("loadDefaultCodonUsage(unknown) produced no frequencies")
	}
}

func TestLoadHumanCodonUsage_Has64Codons(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadHumanCodonUsage()
	// Should have 64 codons (some may be 0 for stop codons)
	if len(model.codonFrequencies) == 0 {
		t.Error("loadHumanCodonUsage() produced no frequencies")
	}
}

func TestLoadHumanCodonUsage_FrequenciesNonNegative(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadHumanCodonUsage()
	for codon, freq := range model.codonFrequencies {
		if freq < 0 {
			t.Errorf("loadHumanCodonUsage: codon %s has negative frequency %f", codon, freq)
		}
	}
}

func TestLoadMouseCodonUsage_HasFrequencies(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadMouseCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadMouseCodonUsage() produced no frequencies")
	}
}

func TestLoadUniversalCodonUsage_HasFrequencies(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadUniversalCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadUniversalCodonUsage() produced no frequencies")
	}
}

// ── Yeast ────────────────────────────────────────────────────────────────────

func TestLoadYeastCodonUsage_HasFrequencies(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadYeastCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadYeastCodonUsage() produced no codon frequencies")
	}
}

func TestLoadYeastCodonUsage_FrequenciesNonNegative(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadYeastCodonUsage()
	for codon, freq := range model.codonFrequencies {
		if freq < 0 {
			t.Errorf("loadYeastCodonUsage: codon %s has negative frequency %f", codon, freq)
		}
	}
}

// Yeast strong codon bias: GAA (Glu) must be more common than GAG.
func TestLoadYeastCodonUsage_GluBias(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadYeastCodonUsage()
	gaa := model.codonFrequencies["GAA"]
	gag := model.codonFrequencies["GAG"]
	if gaa <= gag {
		t.Errorf("yeast Glu bias: GAA (%f) should be > GAG (%f)", gaa, gag)
	}
}

// Yeast: CTG (Leu) must be less common than the universal would give (>25% of Leu
// in universal, but yeast prefers TTA and TTG strongly over CTG).
func TestLoadYeastCodonUsage_LeuBias_CTGNotDominant(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadYeastCodonUsage()
	// In yeast: TTA (26.2‰) > TTG (27.2‰) > CTG (10.0‰)
	// Verify TTA is more frequent than CTG.
	tta := model.codonFrequencies["TTA"]
	ctg := model.codonFrequencies["CTG"]
	if tta <= ctg {
		t.Errorf("yeast Leu bias: TTA (%f) should exceed CTG (%f)", tta, ctg)
	}
}

// Yeast is NOT the same as the universal equal-frequency table.
func TestLoadYeastCodonUsage_DiffersFromUniversal(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadYeastCodonUsage()
	yeastCTG := model.codonFrequencies["CTG"]

	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadUniversalCodonUsage()
	univCTG := model.codonFrequencies["CTG"]

	if yeastCTG == univCTG {
		t.Errorf("yeast CTG (%f) should differ from universal CTG (%f)", yeastCTG, univCTG)
	}
}

// ── E. coli ───────────────────────────────────────────────────────────────────

func TestLoadEcoliCodonUsage_HasFrequencies(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadEcoliCodonUsage()
	if len(model.codonFrequencies) == 0 {
		t.Error("loadEcoliCodonUsage() produced no codon frequencies")
	}
}

func TestLoadEcoliCodonUsage_FrequenciesNonNegative(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadEcoliCodonUsage()
	for codon, freq := range model.codonFrequencies {
		if freq < 0 {
			t.Errorf("loadEcoliCodonUsage: codon %s has negative frequency %f", codon, freq)
		}
	}
}

// E. coli: CTG must dominate Leu (52.4‰ >> all other Leu codons).
func TestLoadEcoliCodonUsage_LeuCTGDominant(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadEcoliCodonUsage()
	ctg := model.codonFrequencies["CTG"]
	tta := model.codonFrequencies["TTA"]
	ttg := model.codonFrequencies["TTG"]
	if ctg <= tta {
		t.Errorf("E. coli Leu bias: CTG (%f) should dominate TTA (%f)", ctg, tta)
	}
	if ctg <= ttg {
		t.Errorf("E. coli Leu bias: CTG (%f) should dominate TTG (%f)", ctg, ttg)
	}
}

// E. coli: GAA (Glu) must be more common than GAG.
func TestLoadEcoliCodonUsage_GluBias(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadEcoliCodonUsage()
	gaa := model.codonFrequencies["GAA"]
	gag := model.codonFrequencies["GAG"]
	if gaa <= gag {
		t.Errorf("E. coli Glu bias: GAA (%f) should be > GAG (%f)", gaa, gag)
	}
}

// E. coli is NOT the same as the universal equal-frequency table.
func TestLoadEcoliCodonUsage_DiffersFromUniversal(t *testing.T) {
	model, _ := NewCodonUsageModel()
	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadEcoliCodonUsage()
	ecoliCTG := model.codonFrequencies["CTG"]

	model.codonFrequencies = make(map[string]float64)
	model.aminoAcidMapping = make(map[string][]string)
	model.loadUniversalCodonUsage()
	univCTG := model.codonFrequencies["CTG"]

	if ecoliCTG == univCTG {
		t.Errorf("E. coli CTG (%f) should differ from universal CTG (%f)", ecoliCTG, univCTG)
	}
}
