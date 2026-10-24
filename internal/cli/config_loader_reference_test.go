//go:build !release

package cli

import (
	"os"
	"path/filepath"
	"testing"
)

// reference_path, reference_region and strict_ref (#27), and variants_path
// and variants_region (#32), load from JSON, flat YAML, v1 YAML and the run
// section of v2 templates.
func TestConfigLoader_ReferenceSettings(t *testing.T) {
	files := map[string]string{
		"c.json": `{"sample_name":"S","read_length":150,"reference_path":"/ref/hg38.fa","reference_region":"chr22:1-1000","strict_ref":true,"variants_path":"calls.vcf.gz","variants_region":"chr5:1-900"}`,
		"flat.yaml": `sample_name: S
read_length: 150
reference_path: /ref/hg38.fa
reference_region: chr22:1-1000
strict_ref: true
variants_path: calls.vcf.gz
variants_region: chr5:1-900
`,
		"v1.yaml": `schema_version: "1.0"
config_type: sequencer
sample:
  name: S
reads:
  length: 150
reference_path: /ref/hg38.fa
reference_region: chr22:1-1000
strict_ref: true
variants_path: calls.vcf.gz
variants_region: chr5:1-900
`,
		"v2.yaml": `version: "2.0"
platform: ont
vendor: Oxford Nanopore Technologies
model: PromethION
specifications:
  read_length: {median: 8000, n50_target: 20000, distribution: log_normal}
run:
  reference_path: /ref/hg38.fa
  reference_region: chr22:1-1000
  strict_ref: true
  variants_path: calls.vcf.gz
  variants_region: chr5:1-900
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
		if c.ReferencePath != "/ref/hg38.fa" || c.ReferenceRegion != "chr22:1-1000" || !c.StrictRef {
			t.Errorf("%s: reference_path %q reference_region %q strict_ref %v", name, c.ReferencePath, c.ReferenceRegion, c.StrictRef)
		}
		if c.VariantsPath != "calls.vcf.gz" || c.VariantsRegion != "chr5:1-900" {
			t.Errorf("%s: variants_path %q variants_region %q", name, c.VariantsPath, c.VariantsRegion)
		}
	}
}
