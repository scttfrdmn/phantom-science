//go:build !release
// +build !release

package cmd

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/spf13/pflag"
)

// setGenerateFlags sets generate flags for one test and restores them after.
func setGenerateFlags(t *testing.T, kv map[string]string) {
	t.Helper()
	t.Cleanup(func() {
		generateCmd.Flags().VisitAll(func(f *pflag.Flag) {
			if f.Changed {
				_ = f.Value.Set(f.DefValue)
				f.Changed = false
			}
		})
	})
	for k, v := range kv {
		if err := generateCmd.Flags().Set(k, v); err != nil {
			t.Fatalf("set --%s: %v", k, err)
		}
	}
}

func writeGenerateConfig(t *testing.T, readLength, reads int) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "seq.json")
	cfg := fmt.Sprintf(`{"sample_name":"CLI","lane":1,"output_dir":%q,"read_length":%d,"reads_per_batch":%d,"total_batches":1,"base_quality":35,"quality_decay":0.1}`,
		filepath.Join(dir, "out"), readLength, reads)
	if err := os.WriteFile(path, []byte(cfg), 0o600); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestGenerateCommand_NewFlags(t *testing.T) {
	for _, name := range []string{"seed", "platform", "read-length-dist", "read-length-n50", "read-length-median",
		"read-length-mean", "read-length-sd", "basecaller-tier", "compress", "region", "total-bases", "total-reads", "workers"} {
		if generateCmd.Flags().Lookup(name) == nil {
			t.Errorf("generateCmd should have --%s flag", name)
		}
	}
}

func TestParseRegion(t *testing.T) {
	r, err := parseRegion("chr7:1,001-2,000")
	if err != nil || r.Chromosome != "chr7" || r.Start != 1000 || r.End != 2000 {
		t.Errorf("parseRegion = %+v, %v", r, err)
	}
	for _, bad := range []string{"chr7", "chr7:10", "chr7:0-5", "chr7:9-5", ":1-5", "chr7:a-b"} {
		if _, err := parseRegion(bad); err == nil {
			t.Errorf("parseRegion(%q) should fail", bad)
		}
	}
}

func TestParseBaseCount(t *testing.T) {
	for in, want := range map[string]int64{"140G": 140e9, "1.5m": 1.5e6, "1000": 1000, "2T": 2e12, "3k": 3000} {
		if got, err := parseBaseCount(in); err != nil || got != want {
			t.Errorf("parseBaseCount(%q) = %d, %v; want %d", in, got, err, want)
		}
	}
	for _, bad := range []string{"", "G", "-5", "x1"} {
		if _, err := parseBaseCount(bad); err == nil {
			t.Errorf("parseBaseCount(%q) should fail", bad)
		}
	}
}

func TestResolveStreamRegion(t *testing.T) {
	vs := []genomics.VariantRecord{
		{Chromosome: "chr1", Position: 500, Ref: "A", Alt: "T"},
		{Chromosome: "chr1", Position: 700, Ref: "ACG", Alt: "A"},
		{Ref: "C", Alt: "G"}, // position-blind: ignored for the span
	}
	r, err := resolveStreamRegion(vs, "", 100)
	if err != nil || r == nil || r.Chromosome != "chr1" || r.Start != 399 || r.End != 802 {
		t.Errorf("span region = %+v, %v; want chr1 [399, 802)", r, err)
	}
	r, _ = resolveStreamRegion(vs[:1], "", 1000)
	if r.Start != 0 {
		t.Errorf("region start not clamped at 0: %+v", r)
	}
	if r, err := resolveStreamRegion(vs[2:], "", 100); r != nil || err != nil {
		t.Errorf("position-blind variants: got %+v, %v; want nil", r, err)
	}
	multi := append(append([]genomics.VariantRecord(nil), vs...), genomics.VariantRecord{Chromosome: "chr2", Position: 5, Ref: "A", Alt: "C"})
	if _, err := resolveStreamRegion(multi, "", 100); err == nil || !strings.Contains(err.Error(), "--region") {
		t.Errorf("multiple chromosomes: err = %v, want a --region hint", err)
	}
	if r, err := resolveStreamRegion(multi, "chr2:1-500", 100); err != nil || r.Chromosome != "chr2" || r.Start != 0 || r.End != 500 {
		t.Errorf("--region: %+v, %v", r, err)
	}
	if _, err := resolveStreamRegion(nil, "chr2:1-50", 100); err == nil {
		t.Error("--region shorter than a read should fail")
	}
}

// generate --stdout with a positioned VCF places reads around the variants
// and writes each allele at its exact coordinate.
func TestGenerateStdout_PositionedVCF(t *testing.T) {
	const readLength = 50
	vcf := filepath.Join(t.TempDir(), "pos.vcf")
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr3\t1000\t.\tA\tAGGGGGGGG\t.\t.\tAF=1\n" +
		"chr3\t1030\t.\tC\tT\t.\t.\tAF=0\n"
	if err := os.WriteFile(vcf, []byte(body), 0o600); err != nil {
		t.Fatal(err)
	}
	setGenerateFlags(t, map[string]string{
		"config": writeGenerateConfig(t, readLength, 400), "stdout": "true", "variants": vcf, "seed": "11",
	})
	var out bytes.Buffer
	if err := runGenerate(generateCmd, &out); err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSpace(out.String()), "\n")
	if len(lines) != 400*4 {
		t.Fatalf("%d lines, want %d: stdout must hold only FASTQ", len(lines), 400*4)
	}
	covered := 0
	for i := 0; i < len(lines); i += 4 {
		header, seq := lines[i], lines[i+1]
		if !strings.HasPrefix(header, "@CLI:1:HXXXXXX:1:1101:") || len(seq) != readLength || strings.Contains(seq, "N") {
			t.Fatalf("bad record %q / %q", header, seq)
		}
		f := strings.Split(strings.Fields(header)[0], ":")
		x, _ := strconv.ParseInt(f[5], 10, 64)
		y, _ := strconv.ParseInt(f[6], 10, 64)
		start := y*100000 + x
		if start < 999-readLength || start >= 1030+readLength {
			t.Fatalf("read start %d outside the variant span", start)
		}
		if off := 999 - start; off >= 0 && off+9 <= readLength {
			covered++
			if seq[off:off+9] != "AGGGGGGGG" {
				t.Errorf("read at %d: insertion missing at offset %d: %s", start, off, seq)
			}
		}
		// The AF=0 SNV at 1030 is always REF ("C"); it sits 8 bases further
		// right in reads that also contain the insertion.
		if off := 1029 - start; off >= 0 {
			if 999-start >= 0 {
				off += 8
			}
			if off < readLength && seq[off] != 'C' {
				t.Errorf("read at %d: REF not enforced at offset %d: %s", start, off, seq)
			}
		}
	}
	if covered < 50 {
		t.Errorf("only %d reads cover the insertion", covered)
	}
}

func TestGenerateStdout_SeedReproducible(t *testing.T) {
	cfg := writeGenerateConfig(t, 60, 50)
	run := func(seed string) []byte {
		setGenerateFlags(t, map[string]string{"config": cfg, "stdout": "true", "seed": seed})
		var out bytes.Buffer
		if err := runGenerate(generateCmd, &out); err != nil {
			t.Fatal(err)
		}
		return out.Bytes()
	}
	if a, b := run("42"), run("42"); !bytes.Equal(a, b) {
		t.Error("--seed 42 runs differ")
	}
}

func TestGenerate_NotImplementedFlags(t *testing.T) {
	for flag, val := range map[string]string{
		"platform": "ont",
		"compress": "gzip", "total-bases": "140G", "total-reads": "1000", "workers": "4",
	} {
		t.Run(flag, func(t *testing.T) {
			setGenerateFlags(t, map[string]string{"config": writeGenerateConfig(t, 50, 5), "stdout": "true", flag: val})
			err := runGenerate(generateCmd, io.Discard)
			if !errors.Is(err, genomics.ErrNotImplemented) {
				t.Errorf("--%s %s: err = %v, want ErrNotImplemented", flag, val, err)
			}
		})
	}
}

func TestGenerate_StdoutOnlyFlags(t *testing.T) {
	for _, flag := range []string{"variants", "region"} {
		setGenerateFlags(t, map[string]string{"config": writeGenerateConfig(t, 50, 5), flag: "x"})
		if err := runGenerate(generateCmd, io.Discard); err == nil || !strings.Contains(err.Error(), "--stdout") {
			t.Errorf("--%s without --stdout: err = %v", flag, err)
		}
	}
}

func TestGenerateCommand(t *testing.T) {
	if generateCmd == nil {
		t.Fatal("generateCmd should not be nil")
	}
	if generateCmd.Use != "generate" {
		t.Errorf("expected Use 'generate', got %q", generateCmd.Use)
	}
	if generateCmd.RunE == nil {
		t.Error("generateCmd should have a RunE function")
	}
	for _, name := range []string{"stdout", "variants", "sample-id", "config"} {
		if generateCmd.Flags().Lookup(name) == nil {
			t.Errorf("generateCmd should have --%s flag", name)
		}
	}
}

func TestGenerateCommand_InheritsModelFlags(t *testing.T) {
	for _, name := range []string{"hmm-model-path", "population", "model-dir"} {
		// Persistent flags from rootCmd appear in InheritedFlags() on child commands.
		if generateCmd.InheritedFlags().Lookup(name) == nil {
			t.Errorf("generateCmd should inherit --%s from rootCmd persistent flags", name)
		}
	}
}

func TestGenerateCommand_MissingConfig(t *testing.T) {
	// Reset the flag to empty before testing.
	if err := generateCmd.Flags().Set("config", ""); err != nil {
		t.Fatalf("reset config flag: %v", err)
	}
	err := generateCmd.RunE(generateCmd, nil)
	if err == nil {
		t.Error("expected error when --config is empty")
	}
	if !strings.Contains(err.Error(), "--config") {
		t.Errorf("error should mention --config, got: %v", err)
	}
}

// TestParseAFFromInfo verifies AF extraction from VCF INFO fields.
func TestParseAFFromInfo(t *testing.T) {
	tests := []struct {
		info    string
		defFreq float64
		want    float64
	}{
		{"AF=0.35;DP=100", 0.5, 0.35},
		{"DP=100;AF=0.1", 0.5, 0.1},
		{".", 0.5, 0.5},        // no AF field → default
		{"AF=1.5", 0.5, 0.5},  // out-of-range → default
		{"AF=0.3,0.2", 0.5, 0.3}, // multi-allelic: use first value
		{"AF=bad", 0.5, 0.5},  // unparseable → default
	}
	for _, tc := range tests {
		got := parseAFFromInfo(tc.info, tc.defFreq)
		if got != tc.want {
			t.Errorf("parseAFFromInfo(%q) = %v, want %v", tc.info, got, tc.want)
		}
	}
}

// TestLoadVariantsVCF verifies variant loading from a synthetic VCF.
func TestLoadVariantsVCF(t *testing.T) {
	vcf := makeSyntheticVCF(t, 5)
	variants, err := loadVariantsVCF(vcf)
	if err != nil {
		t.Fatalf("loadVariantsVCF: %v", err)
	}
	if len(variants) != 5 {
		t.Errorf("expected 5 variants, got %d", len(variants))
	}
	for _, v := range variants {
		if v.Ref == "" || v.Alt == "" {
			t.Errorf("variant has empty Ref or Alt: %+v", v)
		}
		if v.Freq < 0 || v.Freq > 1 {
			t.Errorf("variant Freq %v out of [0,1]", v.Freq)
		}
	}
}

// TestLoadVariantsVCF_Empty verifies an empty VCF returns zero records without error.
func TestLoadVariantsVCF_Empty(t *testing.T) {
	vcf := makeSyntheticVCF(t, 0)
	variants, err := loadVariantsVCF(vcf)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(variants) != 0 {
		t.Errorf("expected 0 variants, got %d", len(variants))
	}
}

// TestLoadVariantsVCF_WithAF verifies that AF= in INFO is parsed correctly.
func TestLoadVariantsVCF_WithAF(t *testing.T) {
	vcfPath := makeTempVCFWithAF(t)
	variants, err := loadVariantsVCF(vcfPath)
	if err != nil {
		t.Fatalf("loadVariantsVCF: %v", err)
	}
	if len(variants) != 2 {
		t.Fatalf("expected 2 variants, got %d", len(variants))
	}
	if variants[0].Freq != 0.35 {
		t.Errorf("variant[0].Freq = %v, want 0.35", variants[0].Freq)
	}
	if variants[1].Freq != 0.1 {
		t.Errorf("variant[1].Freq = %v, want 0.1", variants[1].Freq)
	}
}

func makeTempVCFWithAF(t *testing.T) string {
	t.Helper()
	f, err := os.CreateTemp(t.TempDir(), "af_test_*.vcf")
	if err != nil {
		t.Fatalf("create temp: %v", err)
	}
	fmt.Fprintln(f, "##fileformat=VCFv4.1")
	fmt.Fprintln(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	fmt.Fprintln(f, "chr1\t100\t.\tA\tG\t.\t.\tAF=0.35")
	fmt.Fprintln(f, "chr1\t200\t.\tC\tT\t.\t.\tAF=0.1")
	f.Close()
	return f.Name()
}

// File output is validated when the phantom config (-p) enables validation;
// previously only the root command did this and `generate` ignored it.
func TestGenerate_ValidatesWithPhantomConfig(t *testing.T) {
	writePhantom := func(t *testing.T, validation string) string {
		t.Helper()
		p := filepath.Join(t.TempDir(), "phantom.json")
		if err := os.WriteFile(p, []byte(`{"enterprise":{"customer_id":"C"},"validation":`+validation+`}`), 0o600); err != nil {
			t.Fatal(err)
		}
		return p
	}
	setPhantom := func(t *testing.T, path string) {
		t.Helper()
		if err := rootCmd.PersistentFlags().Set("phantom", path); err != nil {
			t.Fatal(err)
		}
		t.Cleanup(func() { _ = rootCmd.PersistentFlags().Set("phantom", "") })
	}

	setGenerateFlags(t, map[string]string{"config": writeGenerateConfig(t, 100, 200), "seed": "3"})
	setPhantom(t, writePhantom(t, `{"enable":true}`))
	var out bytes.Buffer
	if err := runGenerate(generateCmd, &out); err != nil {
		t.Fatalf("generate with default validation suite: %v\n%s", err, out.String())
	}
	if !strings.Contains(out.String(), "All 9 validation tests passed") {
		t.Errorf("validation did not run on the generated file:\n%s", out.String())
	}

	setGenerateFlags(t, map[string]string{"config": writeGenerateConfig(t, 100, 200), "seed": "3"})
	setPhantom(t, writePhantom(t, `{"enable":true,"statistical_tests":["no_such_test"]}`))
	out.Reset()
	if err := runGenerate(generateCmd, &out); err == nil || !strings.Contains(err.Error(), "failed validation") {
		t.Errorf("failing validation not reported: err=%v\n%s", err, out.String())
	}
}
