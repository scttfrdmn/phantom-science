//go:build !release

package cmd

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/spf13/pflag"
)

// scenario --reference with a paired-end spec: R1 and R2 files of the same
// pairs, every read (errors off) taken from the reference or from the
// reference carrying the planted variant, binned qualities, and a manifest
// that records the reference as the read source.
func TestScenario_ReferencePairedEnd(t *testing.T) {
	ref, seqs := writeCmdReference(t)
	chr1 := seqs["chr1"]
	refBase := chr1[4000]
	alt := "ACGT"[(strings.IndexByte("ACGT", refBase)+1)%4]
	altHap := chr1[:4000] + string(alt) + chr1[4001:]

	dir := t.TempDir()
	r1, r2 := filepath.Join(dir, "p_R1.fastq"), filepath.Join(dir, "p_R2.fastq")
	truth := filepath.Join(dir, "truth.json")
	spec := filepath.Join(dir, "spec.yaml")
	if err := os.WriteFile(spec, []byte(fmt.Sprintf(`phenotype:
  planted_variant: {chrom: chr1, pos: 4001, ref: %c, alt: %c, zygosity: heterozygous}
sequencer: {platform: NovaSeq6000, coverage: 20, read_length: 100, paired_end: true}
output: {proband_fastq: %s, proband_fastq_r2: %s, ground_truth: %s}
`, refBase, alt, r1, r2, truth)), 0o600); err != nil {
		t.Fatal(err)
	}
	cfg := filepath.Join(dir, "seq.json")
	if err := os.WriteFile(cfg, []byte(`{"sample_name":"P","lane":1,"read_length":100,"seed":8,"illumina":{"errors":"none"}}`), 0o600); err != nil {
		t.Fatal(err)
	}

	t.Cleanup(func() {
		for _, fs := range []*pflag.FlagSet{scenarioCmd.Flags(), rootCmd.PersistentFlags()} {
			fs.VisitAll(func(f *pflag.Flag) {
				if f.Changed {
					_ = f.Value.Set(f.DefValue)
					f.Changed = false
				}
			})
		}
	})
	for k, v := range map[string]string{"spec": spec, "reference": ref} {
		if err := scenarioCmd.Flags().Set(k, v); err != nil {
			t.Fatal(err)
		}
	}
	if err := rootCmd.PersistentFlags().Set("config", cfg); err != nil {
		t.Fatal(err)
	}
	if err := scenarioCmd.RunE(scenarioCmd, nil); err != nil {
		t.Fatal(err)
	}

	read := func(path string) (names, seqs []string) {
		f, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		sc := bufio.NewScanner(f)
		for i := 0; sc.Scan(); i++ {
			switch i % 4 {
			case 0:
				names = append(names, strings.Fields(sc.Text())[0])
			case 1:
				seqs = append(seqs, sc.Text())
			case 3:
				for _, q := range sc.Text() {
					if !strings.ContainsRune("#-8F", q) {
						t.Fatalf("%s: quality %q is not a binned Illumina quality", path, q)
					}
				}
			}
		}
		return names, seqs
	}
	n1, s1 := read(r1)
	n2, s2 := read(r2)
	if want := 20 * 8000 / 200; len(s1) != want || len(s2) != want {
		t.Fatalf("got %d R1 and %d R2 reads, want %d pairs", len(s1), len(s2), want)
	}
	for i := range s1 {
		if n1[i] != n2[i] {
			t.Fatalf("pair %d: names %s and %s", i, n1[i], n2[i])
		}
		for _, s := range []string{s1[i], s2[i]} {
			if !strings.Contains(chr1, s) && !strings.Contains(chr1, revcomp(s)) &&
				!strings.Contains(altHap, s) && !strings.Contains(altHap, revcomp(s)) {
				t.Fatalf("read %s is not from the reference: %s", n1[i], s)
			}
		}
	}

	data, err := os.ReadFile(truth)
	if err != nil {
		t.Fatal(err)
	}
	var m struct {
		Reference string `json:"reference"`
		Sequencer struct {
			PairedEnd bool   `json:"paired_end"`
			Bases     string `json:"bases"`
		} `json:"sequencer"`
		ComputeMetrics struct {
			ReadsGenerated int `json:"reads_generated"`
		} `json:"compute_metrics"`
	}
	if err := json.Unmarshal(data, &m); err != nil {
		t.Fatal(err)
	}
	if m.Reference != ref || !m.Sequencer.PairedEnd || m.Sequencer.Bases != "reference" || m.ComputeMetrics.ReadsGenerated != 2*len(s1) {
		t.Fatalf("manifest %+v", m)
	}
}
