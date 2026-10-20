//go:build !release

package clinical

import (
	"encoding/json"
	"fmt"
	"os"
	"time"
)

// GroundTruthManifest is a machine-readable record of what was planted.
type GroundTruthManifest struct {
	GeneratedAt     time.Time         `json:"generated_at"`
	ProbandID       string            `json:"proband_id"`
	PlantedVariants []PlantedVariant  `json:"planted_variants"`
	Trio            *TrioManifest     `json:"trio,omitempty"`
	Sequencer       SequencerManifest `json:"sequencer"`
	ComputeMetrics  ComputeMetrics    `json:"compute_metrics"`
}

// PlantedVariant records a single planted variant in the manifest.
type PlantedVariant struct {
	Gene        string   `json:"gene"`
	Chrom       string   `json:"chrom"`
	Pos         int64    `json:"pos"`
	Ref         string   `json:"ref"`
	Alt         string   `json:"alt"`
	Zygosity    string   `json:"zygosity"`
	ExpectedVAF float64  `json:"expected_vaf"` // 0.5 for het, 1.0 for hom
	Consequence string   `json:"consequence,omitempty"`
	ACMGClass   string   `json:"acmg_class,omitempty"`
	HPOTerms    []string `json:"hpo_terms,omitempty"`
	// Build, ClinVarID and HGVS are set when the variant came from the
	// built-in HPO catalogue.
	Build     string `json:"build,omitempty"`
	ClinVarID string `json:"clinvar_id,omitempty"`
	HGVS      string `json:"hgvs,omitempty"`
}

// TrioManifest records trio generation details.
type TrioManifest struct {
	MotherCarries bool `json:"mother_carries"`
	FatherCarries bool `json:"father_carries"`
	DeNovo        bool `json:"de_novo"`
}

// SequencerManifest records sequencing parameters.
type SequencerManifest struct {
	Platform      string `json:"platform"`
	Coverage      int    `json:"coverage"`
	ReadLength    int    `json:"read_length"`
	RegionFocused string `json:"region_focused,omitempty"` // e.g. "chr1:42904241-42954241"
}

// ComputeMetrics records generation performance metrics.
type ComputeMetrics struct {
	ReadsGenerated int   `json:"reads_generated"`
	WallClockNS    int64 `json:"wall_clock_ns"`
}

// WriteManifest serialises the manifest as indented JSON and writes it to path.
func WriteManifest(path string, m *GroundTruthManifest) error {
	data, err := json.MarshalIndent(m, "", "  ")
	if err != nil {
		return fmt.Errorf("marshalling manifest: %w", err)
	}
	if err := os.WriteFile(path, data, 0o600); err != nil { // #nosec G306
		return fmt.Errorf("writing manifest to %s: %w", path, err)
	}
	return nil
}
