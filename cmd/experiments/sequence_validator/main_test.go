//go:build !release

package main

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestParseFilename(t *testing.T) {
	cases := []struct {
		name    string
		pop     string
		length  int
		seq     int
		wantErr bool
	}{
		{"YRI_L1000_S01.fasta", "YRI", 1000, 1, false},
		{"YRI_L50_S7.fa", "YRI", 50, 7, false},
		{"HG_YRI_L200_S03.FASTA", "HG_YRI", 200, 3, false}, // underscore in population
		{"invalid_name.fasta", "", 0, 0, true},
		{"YRI_1000_S01.fasta", "", 0, 0, true}, // missing L
		{"YRI_L1000_01.fasta", "", 0, 0, true}, // missing S
		{"YRI_Lxxx_S01.fasta", "", 0, 0, true},
		{"YRI_L0_S01.fasta", "", 0, 0, true},
		{"YRI_L1000_Sxx.fasta", "", 0, 0, true},
		{"_L1000_S01.fasta", "", 0, 0, true},
	}
	for _, c := range cases {
		pop, length, seq, err := parseFilename(c.name)
		if c.wantErr {
			if err == nil {
				t.Errorf("parseFilename(%q) = %q,%d,%d; want error", c.name, pop, length, seq)
			}
			continue
		}
		if err != nil || pop != c.pop || length != c.length || seq != c.seq {
			t.Errorf("parseFilename(%q) = %q,%d,%d,%v; want %q,%d,%d", c.name, pop, length, seq, err, c.pop, c.length, c.seq)
		}
	}
}

func TestParseHeaderMetadata(t *testing.T) {
	gc, ok, region := parseHeaderMetadata("YRI_L1000_S01 population=YRI region=Africa model_target_gc=0.496")
	if !ok || math.Abs(gc-0.496) > 1e-12 || region != "Africa" {
		t.Errorf("got gc=%v ok=%t region=%q", gc, ok, region)
	}
	if _, ok, _ := parseHeaderMetadata("YRI_L1000_S01 population=YRI"); ok {
		t.Error("hasTarget true for a header without model_target_gc")
	}
	if _, ok, _ := parseHeaderMetadata("x model_target_gc=abc"); ok {
		t.Error("hasTarget true for a non-numeric model_target_gc")
	}
}

func TestCountValidBases(t *testing.T) {
	if v, i := countValidBases("ATCGatcgNNXX"); v != 8 || i != 4 {
		t.Errorf("countValidBases = %d,%d; want 8,4", v, i)
	}
}

func TestCalculateGCContent(t *testing.T) {
	cases := map[string]float64{
		"GCGCGCGC": 1,
		"ATATATAT": 0,
		"GCATgcat": 0.5,
		"":         0,
		"NNNN":     0,
		// N is not a base of known composition: 2 GC of 4 ACGT bases, not 2 of 6.
		"GGNNAA": 0.5,
	}
	for seq, want := range cases {
		if got := calculateGCContent(seq); math.Abs(got-want) > 1e-12 {
			t.Errorf("calculateGCContent(%q) = %v, want %v", seq, got, want)
		}
	}
}

func TestCalculateComplexity(t *testing.T) {
	// 7 windows of 4, distinct: ACGT CGTA GTAC TACG.
	if got := calculateComplexity("ACGTACGTAC", 4); math.Abs(got-4.0/7) > 1e-12 {
		t.Errorf("complexity = %v, want 4/7", got)
	}
	if got := calculateComplexity(strings.Repeat("A", 100), 4); math.Abs(got-1.0/97) > 1e-12 {
		t.Errorf("homopolymer complexity = %v, want 1/97", got)
	}
	// Windows spanning N are not k-mers: ACGT and TTTT only.
	if got := calculateComplexity("ACGTNTTTT", 4); got != 1 {
		t.Errorf("complexity with N = %v, want 1 (2 distinct of 2)", got)
	}
	if got := calculateComplexity("acgt", 4); got != 1 {
		t.Errorf("lower-case complexity = %v, want 1", got)
	}
	if got := calculateComplexity("AT", 8); got != 0 {
		t.Errorf("complexity of a sequence shorter than k = %v, want 0", got)
	}
}

func TestFindMaxHomopolymer(t *testing.T) {
	cases := map[string]int{"": 0, "ATCG": 1, "ATAAAATCG": 4, "ATCGGGGG": 5, "aAAa": 4, "GGGGAAAAAAT": 6}
	for seq, want := range cases {
		if got := findMaxHomopolymer(seq); got != want {
			t.Errorf("findMaxHomopolymer(%q) = %d, want %d", seq, got, want)
		}
	}
}

func TestClassifyComplexity(t *testing.T) {
	cases := map[float64]string{0.9: "Very High", 0.8: "Very High", 0.7: "High", 0.5: "Medium", 0.3: "Low", 0.1: "Very Low"}
	for c, want := range cases {
		if got := classifyComplexity(c); got != want {
			t.Errorf("classifyComplexity(%v) = %q, want %q", c, got, want)
		}
	}
}

func TestReadSingleFASTA(t *testing.T) {
	h, s, err := readSingleFASTA(strings.NewReader(">id a=b\nACGT\n  acgt  \n\n"))
	if err != nil || h != "id a=b" || s != "ACGTacgt" {
		t.Errorf("got %q %q %v", h, s, err)
	}
	for name, in := range map[string]string{
		"empty":          "",
		"no header":      "ACGT\n",
		"empty sequence": ">id\n\n",
		"two records":    ">a\nACGT\n>b\nACGT\n",
	} {
		if _, _, err := readSingleFASTA(strings.NewReader(in)); err == nil {
			t.Errorf("%s: no error", name)
		}
	}
}

// bufio.Scanner's default 64 KiB token limit used to reject sequences
// written on one line.
func TestReadSingleFASTA_LongLine(t *testing.T) {
	seq := strings.Repeat("ACGT", 100000)
	_, s, err := readSingleFASTA(strings.NewReader(">id\n" + seq + "\n"))
	if err != nil || len(s) != len(seq) {
		t.Fatalf("len=%d err=%v; want %d bases", len(s), err, len(seq))
	}
}

// randomSeq returns seeded uniform ACGT sequence and its GC fraction.
func randomSeq(n int, seed uint64) (string, float64) {
	rng := rand.New(rand.NewPCG(seed, 1))
	b := make([]byte, n)
	gc := 0
	for i := range b {
		b[i] = "ACGT"[rng.IntN(4)]
		if b[i] == 'C' || b[i] == 'G' {
			gc++
		}
	}
	return string(b), float64(gc) / float64(n)
}

func writeFASTA(t *testing.T, dir, name, header, seq string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(">"+header+"\n"+seq+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestValidateSequenceFile_PassingSequence(t *testing.T) {
	seq, gc := randomSeq(1000, 1)
	p := writeFASTA(t, t.TempDir(), "YRI_L1000_S01.fasta",
		fmt.Sprintf("YRI_L1000_S01 region=Africa model_target_gc=%.4f", gc), seq)
	r, err := validateFile(p, defaultThresholds)
	if err != nil {
		t.Fatal(err)
	}
	if !r.IsValid || len(r.Issues) != 0 {
		t.Errorf("random 1 kb sequence at its own GC failed: %v", r.Issues)
	}
	if r.Population != "YRI" || r.Region != "Africa" || r.ActualLength != 1000 || !r.HasTargetGC {
		t.Errorf("unexpected result %+v", r)
	}
}

// Without model_target_gc the target used to be 0, so every sequence was
// flagged for a ~0.5 GC deviation.
func TestValidateSequenceFile_NoTargetGCIsNotAnIssue(t *testing.T) {
	seq, _ := randomSeq(1000, 2)
	p := writeFASTA(t, t.TempDir(), "CEU_L1000_S01.fasta", "CEU_L1000_S01 region=Europe", seq)
	r, err := validateFile(p, defaultThresholds)
	if err != nil {
		t.Fatal(err)
	}
	if r.HasTargetGC || !r.IsValid {
		t.Errorf("HasTargetGC=%t IsValid=%t issues=%v; want no target and valid", r.HasTargetGC, r.IsValid, r.Issues)
	}
}

func TestValidateSequenceFile_Failures(t *testing.T) {
	dir := t.TempDir()
	seq, gc := randomSeq(1000, 3)
	cases := []struct {
		name, header, seq, issue string
	}{
		{"A_L1000_S01.fasta", fmt.Sprintf("x model_target_gc=%.3f", gc+0.1), seq, "differs from target"},
		{"A_L999_S02.fasta", "x", seq, "length 1000, file name says 999"},
		{"A_L1000_S03.fasta", "x", seq[:500] + "NN" + seq[502:], "2 non-ACGT bases"},
		{"A_L1000_S04.fasta", "x", seq[:400] + strings.Repeat("G", 12) + seq[412:], "homopolymer run of 12"},
		{"A_L1000_S05.fasta", "x", strings.Repeat("ACGTTGCA", 125), "complexity"},
	}
	for _, c := range cases {
		r, err := validateFile(writeFASTA(t, dir, c.name, c.header, c.seq), defaultThresholds)
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		if r.IsValid || !strings.Contains(strings.Join(r.Issues, "; "), c.issue) {
			t.Errorf("%s: IsValid=%t issues=%v; want an issue containing %q", c.name, r.IsValid, r.Issues, c.issue)
		}
	}
}

func TestCalculateValidationScore(t *testing.T) {
	perfect := &ValidationResult{Complexity: 0.8}
	if s := calculateValidationScore(perfect); s != 1 {
		t.Errorf("perfect score = %v, want 1 (clamped)", s)
	}
	// 1 - 2*0.1 = 0.8; the deviation only counts when there is a target.
	dev := &ValidationResult{Complexity: 0.5, HasTargetGC: true, GCDeviation: -0.1}
	if s := calculateValidationScore(dev); math.Abs(s-0.8) > 1e-12 {
		t.Errorf("score with GC deviation = %v, want 0.8", s)
	}
	dev.HasTargetGC = false
	if s := calculateValidationScore(dev); s != 1 {
		t.Errorf("score ignoring absent target = %v, want 1", s)
	}
	worst := &ValidationResult{Complexity: 0.1, InvalidBases: 1, HasTargetGC: true, GCDeviation: 0.5,
		Issues: []string{"a", "b", "c"}}
	if s := calculateValidationScore(worst); s != 0 {
		t.Errorf("worst score = %v, want 0 (clamped)", s)
	}
}

func TestGenerateValidationReport_Aggregates(t *testing.T) {
	results := []ValidationResult{
		{Population: "YRI", Region: "Africa", IsValid: true, ValidationScore: 0.9, GCContent: 0.40, Complexity: 0.8},
		{Population: "CEU", Region: "Europe", IsValid: false, ValidationScore: 0.4, GCContent: 0.50, Complexity: 0.6, Issues: []string{"x"}},
		{Population: "YRI", Region: "Africa", IsValid: false, ValidationScore: 0.7, GCContent: 0.44, Complexity: 0.6, Issues: []string{"x"}},
		{Population: "PEL", Region: "Americas", IsValid: true, ValidationScore: 0.4},
	}
	r := generateValidationReport(results, defaultThresholds)
	if r.TotalSequences != 4 || r.ValidSequences != 2 {
		t.Errorf("total/valid = %d/%d, want 4/2", r.TotalSequences, r.ValidSequences)
	}
	if math.Abs(r.OverallScore-0.6) > 1e-12 {
		t.Errorf("overall score = %v, want 0.6", r.OverallScore)
	}
	var order []string
	for _, s := range r.PopulationSummary {
		order = append(order, s.Population)
	}
	// YRI 0.8 first; CEU and PEL tie at 0.4 and sort by name.
	if strings.Join(order, ",") != "YRI,CEU,PEL" {
		t.Errorf("population order %v, want YRI,CEU,PEL", order)
	}
	yri := r.PopulationSummary[0]
	if yri.TotalSequences != 2 || yri.ValidSequences != 1 || math.Abs(yri.AvgGCContent-0.42) > 1e-12 ||
		math.Abs(yri.AvgComplexity-0.7) > 1e-12 || yri.Issues["x"] != 1 {
		t.Errorf("YRI summary %+v", yri)
	}
	if math.Abs(r.RegionalAnalysis["Africa"]-0.8) > 1e-12 || r.RegionalAnalysis["Americas"] != 0.4 {
		t.Errorf("regional analysis %v", r.RegionalAnalysis)
	}
	if r.Timestamp == "" {
		t.Error("timestamp not set")
	}
}

// The display used to list only four hard-coded regions (one of them,
// "East Asia", impossible to parse from a header) and divide by zero on an
// empty report.
func TestDisplayValidationResults(t *testing.T) {
	var out bytes.Buffer
	displayValidationResults(&out, &ValidationReport{
		TotalSequences: 1, ValidSequences: 1,
		PopulationSummary: []PopulationSummary{{Population: "PNG", Region: "Oceania", TotalSequences: 1, ValidSequences: 1,
			Issues: map[string]int{"b issue": 2, "a issue": 2, "c issue": 5}}},
		RegionalAnalysis: map[string]float64{"Oceania": 0.9, "EastAsia": 0.8},
	})
	s := out.String()
	for _, want := range []string{"Oceania: 0.900", "EastAsia: 0.800", "PNG (Oceania)"} {
		if !strings.Contains(s, want) {
			t.Errorf("output missing %q:\n%s", want, s)
		}
	}
	if !(strings.Index(s, "c issue") < strings.Index(s, "a issue") && strings.Index(s, "a issue") < strings.Index(s, "b issue")) {
		t.Errorf("issues not ordered by count then name:\n%s", s)
	}

	out.Reset()
	displayValidationResults(&out, generateValidationReport(nil, defaultThresholds))
	if strings.Contains(out.String(), "NaN") {
		t.Errorf("empty report printed NaN:\n%s", out.String())
	}
}

func TestRun_ExitStatusAndReport(t *testing.T) {
	dir := t.TempDir()
	good, gc := randomSeq(1000, 10)
	writeFASTA(t, dir, "YRI_L1000_S01.fasta", fmt.Sprintf("h region=Africa model_target_gc=%.4f", gc), good)
	sub := filepath.Join(dir, "CEU")
	if err := os.Mkdir(sub, 0o700); err != nil {
		t.Fatal(err)
	}
	writeFASTA(t, sub, "CEU_L1000_S01.fasta", "h region=Europe", good[:900]) // wrong length
	writeFASTA(t, dir, "notes.fasta", "h", good)                             // unparseable name
	if err := os.WriteFile(filepath.Join(dir, "README.txt"), []byte("x"), 0o600); err != nil {
		t.Fatal(err)
	}

	report := filepath.Join(t.TempDir(), "report.json")
	var stdout, stderr bytes.Buffer
	if code := run([]string{"-out", report, dir}, &stdout, &stderr); code != 2 {
		t.Fatalf("exit %d, want 2 (one sequence fails)\nstdout:\n%s\nstderr:\n%s", code, stdout.String(), stderr.String())
	}
	data, err := os.ReadFile(report)
	if err != nil {
		t.Fatal(err)
	}
	var got ValidationReport
	if err := json.Unmarshal(data, &got); err != nil {
		t.Fatal(err)
	}
	if got.TotalSequences != 2 || got.ValidSequences != 1 {
		t.Errorf("report total/valid = %d/%d, want 2/1", got.TotalSequences, got.ValidSequences)
	}
	if len(got.Skipped) != 1 || filepath.Base(got.Skipped[0].Path) != "notes.fasta" {
		t.Errorf("skipped = %+v, want notes.fasta", got.Skipped)
	}
	if !strings.Contains(stderr.String(), "skipped") {
		t.Errorf("stderr does not mention the skipped file:\n%s", stderr.String())
	}

	// Just the passing file, given as a file path.
	if code := run([]string{"-out", report, filepath.Join(dir, "YRI_L1000_S01.fasta")}, &stdout, &stderr); code != 0 {
		t.Errorf("exit %d for a passing file, want 0", code)
	}
	// Tightening a threshold makes it fail.
	if code := run([]string{"-out", report, "-max-homopolymer", "2", filepath.Join(dir, "YRI_L1000_S01.fasta")}, &stdout, &stderr); code != 2 {
		t.Errorf("exit %d with -max-homopolymer 2, want 2", code)
	}
}

func TestRun_UsageErrors(t *testing.T) {
	var stdout, stderr bytes.Buffer
	if code := run(nil, &stdout, &stderr); code != 1 || !strings.Contains(stderr.String(), "usage") {
		t.Errorf("no paths: exit %d, stderr %q; want 1 and usage", code, stderr.String())
	}
	report := filepath.Join(t.TempDir(), "r.json")
	if code := run([]string{"-out", report, t.TempDir()}, &stdout, &stderr); code != 1 {
		t.Errorf("empty directory: exit %d, want 1", code)
	}
	if code := run([]string{"-out", report, filepath.Join(t.TempDir(), "missing")}, &stdout, &stderr); code != 1 {
		t.Errorf("missing path: exit %d, want 1", code)
	}
}

func TestSaveValidationReport_InvalidPath(t *testing.T) {
	if err := saveValidationReport(&ValidationReport{}, filepath.Join(t.TempDir(), "no", "such", "r.json")); err == nil {
		t.Error("expected an error writing into a missing directory")
	}
}
