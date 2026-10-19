//go:build !release

package integration

import (
	"os"
	"path/filepath"
	"testing"
)

// The tool-output parsers must report what the tool printed and fail when
// they find nothing to parse, never substitute a plausible default.

const flagstatSample = `200 + 0 in total (QC-passed reads + QC-failed reads)
200 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
197 + 0 mapped (98.50% : N/A)
197 + 0 primary mapped (98.50% : N/A)
`

func TestParseMappingRate(t *testing.T) {
	got, err := parseMappingRate(flagstatSample)
	if err != nil || got != 0.985 {
		t.Errorf("parseMappingRate = %v, %v; want 0.985, nil", got, err)
	}
	for name, in := range map[string]string{
		"empty":   "",
		"garbage": "samtools: command failed\n",
		"no rate": "200 + 0 in total (QC-passed reads + QC-failed reads)\n",
	} {
		if got, err := parseMappingRate(in); err == nil {
			t.Errorf("%s: parseMappingRate = %v, nil; want an error", name, got)
		}
	}
}

const bcftoolsStatsSample = "# SN\t[2]id\t[3]key\t[4]value\n" +
	"SN\t0\tnumber of samples:\t0\n" +
	"SN\t0\tnumber of records:\t40\n" +
	"SN\t0\tnumber of SNPs:\t30\n" +
	"SN\t0\tnumber of indels:\t10\n" +
	"# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\n" +
	"TSTV\t0\t20\t10\t2.00\t20\t10\t2.00\n"

func TestParseBCFtoolsStats(t *testing.T) {
	st, err := parseBCFtoolsStats(bcftoolsStatsSample)
	if err != nil {
		t.Fatal(err)
	}
	want := VariantStatistics{TotalVariants: 40, SNPCount: 30, IndelCount: 10, TransitionCount: 20, TransversionCount: 10}
	if st != want {
		t.Errorf("parseBCFtoolsStats = %+v, want %+v", st, want)
	}
	for name, in := range map[string]string{
		"empty":   "",
		"no TSTV": "SN\t0\tnumber of records:\t40\n",
		"garbage": "[E::main] failed to open file\n",
	} {
		if st, err := parseBCFtoolsStats(in); err == nil {
			t.Errorf("%s: parseBCFtoolsStats = %+v, nil; want an error", name, st)
		}
	}
}

func TestParseVCFStatistics(t *testing.T) {
	dir := t.TempDir()
	vcf := filepath.Join(dir, "x.vcf")
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\n" +
		"chr1\t200\t.\tAT\tA\t50\tLowQual\t.\n"
	if err := os.WriteFile(vcf, []byte(body), 0o600); err != nil {
		t.Fatal(err)
	}
	st, err := parseVCFStatistics(vcf)
	if err != nil {
		t.Fatal(err)
	}
	want := VCFStatistics{VariantCount: 2, SNPCount: 1, IndelCount: 1, QualityScore: 40, FilterPassRate: 0.5}
	if st != want {
		t.Errorf("parseVCFStatistics = %+v, want %+v", st, want)
	}

	if _, err := parseVCFStatistics(filepath.Join(dir, "missing.vcf")); err == nil {
		t.Error("parseVCFStatistics of a missing file returned no error")
	}
	notVCF := filepath.Join(dir, "not.vcf")
	if err := os.WriteFile(notVCF, []byte("not a vcf\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	if _, err := parseVCFStatistics(notVCF); err == nil {
		t.Error("parseVCFStatistics of a file without a #CHROM header returned no error")
	}
}

func TestParseMeanMAPQ(t *testing.T) {
	sam := "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\nr2\t0\tchr1\t5\t20\t4M\t*\t0\t0\tACGT\tIIII\n"
	if got, err := parseMeanMAPQ(sam); err != nil || got != 40 {
		t.Errorf("parseMeanMAPQ = %v, %v; want 40, nil", got, err)
	}
	if _, err := parseMeanMAPQ(""); err == nil {
		t.Error("parseMeanMAPQ of no records returned no error")
	}
}
