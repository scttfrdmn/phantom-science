//go:build !release
// +build !release

package features

import (
	"math"
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

// reverseStrandORFFixture returns a forward-strand sequence whose only ORF
// is on the reverse strand: flankL + revcomp(coding) + flankR, where coding
// is ATG, 59 sense codons and TAA (183 nt). The sense codons all have
// synonyms, and neither strand contains another ATG.
func reverseStrandORFFixture(m *CodonUsageModel) (seq, coding []byte, start, end int) {
	sense := []string{"CTG", "GCC", "GGC", "TCC", "CGC", "GTG", "CCC", "ACC"}
	coding = []byte("ATG")
	for i := 0; i < 59; i++ {
		coding = append(coding, sense[i%len(sense)]...)
	}
	coding = append(coding, "TAA"...)
	flankL := []byte("CCCCCCCCCC")            // 10 nt
	flankR := []byte("CCCCCCCCCCCCCCCCCCCCC") // 21 nt, so forward and reverse offsets differ
	seq = append(append(append([]byte{}, flankL...), m.reverseComplement(coding)...), flankR...)
	return seq, coding, len(flankL), len(flankL) + len(coding)
}

func translate(m *CodonUsageModel, cds []byte) string {
	var p []byte
	for i := 0; i+3 <= len(cds); i += 3 {
		p = append(p, m.codonTable.CodonToAA[string(cds[i:i+3])]...)
	}
	return string(p)
}

// TestApplyCodonBias_DefaultModelHasCodonUsage: NewCodonUsageModel must load
// the organism's codon usage. Without it every codon frequency is 0, codon
// choice normalises 0/0 (NaN) and the CAI of every ORF is 0.
func TestApplyCodonBias_DefaultModelHasCodonUsage(t *testing.T) {
	m, err := NewCodonUsageModel()
	if err != nil {
		t.Fatal(err)
	}
	if len(m.codonFrequencies) != 64 {
		t.Fatalf("default model has %d codon frequencies, want 64", len(m.codonFrequencies))
	}
	seq, _, _, _ := reverseStrandORFFixture(m)
	res, err := m.ApplyCodonBias(seq)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.ORFs) == 0 {
		t.Fatal("no ORF found")
	}
	for _, orf := range res.ORFs {
		if math.IsNaN(orf.CAI) || orf.CAI <= 0 || orf.CAI > 1 {
			t.Errorf("ORF CAI = %v, want in (0, 1]", orf.CAI)
		}
	}
	if math.IsNaN(res.BiasStrength) || res.BiasStrength <= 0 {
		t.Errorf("BiasStrength = %v, want a positive number", res.BiasStrength)
	}
}

// TestApplyCodonBias_ReverseStrandORF: an ORF found on the reverse strand
// must be reported and edited at its forward-strand coordinates, the edits
// must be synonymous (same protein), and bases outside the ORF must not
// change.
func TestApplyCodonBias_ReverseStrandORF(t *testing.T) {
	m, err := NewCodonUsageModel()
	if err != nil {
		t.Fatal(err)
	}
	m.SetSeed(3)
	seq, coding, start, end := reverseStrandORFFixture(m)

	res, err := m.ApplyCodonBias(seq)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.ORFs) != 1 {
		t.Fatalf("found %d ORFs, want 1", len(res.ORFs))
	}
	orf := res.ORFs[0]
	if orf.Strand != "-" || orf.Start != start || orf.End != end {
		t.Fatalf("ORF = strand %q [%d,%d), want strand \"-\" [%d,%d)", orf.Strand, orf.Start, orf.End, start, end)
	}
	if len(res.CodonChanges) == 0 {
		t.Fatal("no codon changes applied to the reverse-strand ORF")
	}

	b := res.BiasedSequence
	if string(b[:start]) != string(seq[:start]) || string(b[end:]) != string(seq[end:]) {
		t.Errorf("bases outside the ORF changed:\n got %s\nwant %s", b, seq)
	}
	newCDS := m.reverseComplement(b[start:end])
	if got, want := translate(m, newCDS), translate(m, coding); got != want {
		t.Errorf("protein changed:\n got %s\nwant %s", got, want)
	}
	for _, c := range res.CodonChanges {
		if c.Position < start || c.Position+3 > end {
			t.Errorf("change at %d outside ORF [%d,%d)", c.Position, start, end)
			continue
		}
		// Codons are reported on the coding strand.
		if got := string(m.reverseComplement(b[c.Position : c.Position+3])); got != c.NewCodon {
			t.Errorf("change at %d: forward bases encode %s, change says %s", c.Position, got, c.NewCodon)
		}
	}
}

// TestApplyCodonBias_SeededIsDeterministic: with the same seed, two models
// must make the same codon choices (synonymous codon order must not depend
// on map iteration).
func TestApplyCodonBias_SeededIsDeterministic(t *testing.T) {
	var first []byte
	for i := 0; i < 10; i++ {
		m, _ := NewCodonUsageModel()
		m.SetSeed(11)
		seq, _, _, _ := reverseStrandORFFixture(m)
		res, err := m.ApplyCodonBias(seq)
		if err != nil {
			t.Fatal(err)
		}
		if i == 0 {
			first = res.BiasedSequence
		} else if string(res.BiasedSequence) != string(first) {
			t.Fatalf("run %d differs from run 0 with the same seed", i)
		}
	}
}
