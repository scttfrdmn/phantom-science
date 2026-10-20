//go:build !release

package cli

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// The illumina section (#27) loads from JSON, flat YAML and v1 YAML.
func TestConfigLoader_IlluminaSection(t *testing.T) {
	want := types.IlluminaSpec{
		QualityModel: types.IlluminaQualityBinned, Errors: types.IlluminaErrorsNone,
		InsertSize: types.InsertSizeSpec{Dist: types.InsertSizeLognormal, Mean: 450, SD: 90, Min: 200, Max: 900},
	}
	files := map[string]string{
		"c.json": `{"sample_name":"S","read_length":150,"illumina":{"quality_model":"binned","errors":"none",
			"insert_size":{"dist":"lognormal","mean":450,"sd":90,"min":200,"max":900}}}`,
		"flat.yaml": `sample_name: S
read_length: 150
illumina:
  quality_model: binned
  errors: none
  insert_size: {dist: lognormal, mean: 450, sd: 90, min: 200, max: 900}
`,
		"v1.yaml": `schema_version: "1.0"
config_type: sequencer
sample:
  name: S
reads:
  length: 150
illumina:
  quality_model: binned
  errors: none
  insert_size:
    dist: lognormal
    mean: 450
    sd: 90
    min: 200
    max: 900
`,
	}
	for name, content := range files {
		path := filepath.Join(t.TempDir(), name)
		if err := os.WriteFile(path, []byte(content), 0o600); err != nil {
			t.Fatal(err)
		}
		c, err := NewConfigLoader(false).LoadSequencerConfig(path)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if c.Illumina != want {
			t.Errorf("%s: illumina = %+v, want %+v", name, c.Illumina, want)
		}
	}
}
