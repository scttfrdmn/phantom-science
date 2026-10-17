//go:build !release
// +build !release

package cmd

import (
	"fmt"
	"os"
	"strings"
	"testing"
)

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
