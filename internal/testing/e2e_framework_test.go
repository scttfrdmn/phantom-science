//go:build !release
// +build !release

package testing

import (
	"bufio"
	"bytes"
	"errors"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"runtime"
	"strconv"
	"strings"
	"testing"
	"time"
)

// modeledPlatforms are the profiles the engine generates reads for.
var modeledPlatforms = []SequencingPlatform{IlluminaNovaSeq6000, IlluminaNextSeq2000, IlluminaMiSeq, ONTPromethION24}

// generate returns data for platform with the given seed and read count,
// failing the test on error.
func generate(t *testing.T, f *E2EFramework, platform SequencingPlatform, seed uint64, reads int) *SyntheticData {
	t.Helper()
	c := f.CreatePlatformConfig(platform, "CEU")
	if c == nil {
		t.Fatalf("no profile for %s", platform)
	}
	c.Seed = seed
	if reads > 0 {
		c.Reads = reads
	}
	data, err := f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatalf("GenerateSyntheticData(%s): %v", platform, err)
	}
	return data
}

func TestCreatePlatformConfig_Unknown(t *testing.T) {
	if c := NewE2EFramework().CreatePlatformConfig("no_such_platform", "CEU"); c != nil {
		t.Errorf("CreatePlatformConfig(unknown) = %+v, want nil", c)
	}
}

// TestCreatePlatformConfig_NoAliasing is a regression test: the config was
// copied with `config := *baseConfig`, so every returned config shared the
// profile's TechnologyParams map and writing to one changed the framework's
// profile and every other config.
func TestCreatePlatformConfig_NoAliasing(t *testing.T) {
	f := NewE2EFramework()
	a := f.CreatePlatformConfig(IlluminaNovaSeq6000, "YRI")
	b := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	want := b.TechnologyParams["chemistry"]
	if want == "" {
		t.Fatal("NovaSeq profile has no chemistry label")
	}

	a.TechnologyParams["chemistry"] = "mutated"
	a.TechnologyParams["added"] = "x"
	a.ReadLength = 1

	if got := b.TechnologyParams["chemistry"]; got != want {
		t.Errorf("mutating one config changed another: chemistry = %q, want %q", got, want)
	}
	c := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CHB")
	if got := c.TechnologyParams["chemistry"]; got != want {
		t.Errorf("mutating a config changed the framework profile: chemistry = %q, want %q", got, want)
	}
	if _, ok := c.TechnologyParams["added"]; ok {
		t.Error("key added to one config appears in a new config")
	}
	if c.ReadLength != 150 {
		t.Errorf("ReadLength = %d, want 150", c.ReadLength)
	}
	if b.Population != "CEU" || a.Population != "YRI" {
		t.Errorf("populations = %q, %q; want YRI, CEU", a.Population, b.Population)
	}
}

// TestGenerateSyntheticData_NotModeled checks that platforms without a read
// model return ErrPlatformNotModeled instead of fabricated reads.
func TestGenerateSyntheticData_NotModeled(t *testing.T) {
	f := NewE2EFramework()
	for _, p := range []SequencingPlatform{PacBioRevio, MGI_DNBSEQ_T20x2, ElementAVITI} {
		if f.Modeled(p) {
			t.Errorf("%s reported as modeled", p)
		}
		data, err := f.GenerateSyntheticData(f.CreatePlatformConfig(p, "CEU"))
		if !errors.Is(err, ErrPlatformNotModeled) {
			t.Errorf("%s: err = %v, want ErrPlatformNotModeled", p, err)
		}
		if data != nil {
			t.Errorf("%s: returned data with the error", p)
		}
		if _, err := f.RunBenchmark(f.CreatePlatformConfig(p, "CEU")); !errors.Is(err, ErrPlatformNotModeled) {
			t.Errorf("%s: RunBenchmark err = %v, want ErrPlatformNotModeled", p, err)
		}
	}
	for _, p := range modeledPlatforms {
		if !f.Modeled(p) {
			t.Errorf("%s reported as not modeled", p)
		}
	}
}

func TestGenerateSyntheticData_InvalidConfig(t *testing.T) {
	f := NewE2EFramework()
	if _, err := f.GenerateSyntheticData(nil); err == nil {
		t.Error("nil config: want error")
	}
	c := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	c.Seed = 0
	if _, err := f.GenerateSyntheticData(c); err == nil {
		t.Error("zero seed: want error (the run could not be reproduced)")
	}
	c = f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	c.Reads = 0
	if _, err := f.GenerateSyntheticData(c); err == nil {
		t.Error("zero reads: want error")
	}
}

// TestGenerateSyntheticData_MatchesConfig checks the engine's reads against
// the profile: read count, read lengths, quality encoding, and that Stats
// is a recomputation from the reads.
func TestGenerateSyntheticData_MatchesConfig(t *testing.T) {
	f := NewE2EFramework()
	for _, p := range modeledPlatforms {
		t.Run(string(p), func(t *testing.T) {
			c := f.CreatePlatformConfig(p, "YRI")
			data, err := f.GenerateSyntheticData(c)
			if err != nil {
				t.Fatal(err)
			}
			if len(data.Sequences) != c.Reads || len(data.QualityScores) != c.Reads || len(data.Headers) != c.Reads {
				t.Fatalf("got %d/%d/%d sequences/qualities/headers, want %d",
					len(data.Sequences), len(data.QualityScores), len(data.Headers), c.Reads)
			}
			if data.Platform != p || data.Population != "YRI" || data.Seed != c.Seed {
				t.Errorf("labels = %s/%s/%d, want %s/YRI/%d", data.Platform, data.Population, data.Seed, p, c.Seed)
			}
			spec := c.ReadLengthDistribution
			for i, seq := range data.Sequences {
				if spec.Dist == "" {
					if len(seq) != c.ReadLength {
						t.Fatalf("read %d length %d, want %d", i, len(seq), c.ReadLength)
					}
				} else if len(seq) < spec.Min || len(seq) > spec.Max {
					t.Fatalf("read %d length %d outside [%d, %d]", i, len(seq), spec.Min, spec.Max)
				}
				if len(data.QualityScores[i]) != len(seq) {
					t.Fatalf("read %d: %d qualities for %d bases", i, len(data.QualityScores[i]), len(seq))
				}
				for _, q := range data.QualityScores[i] {
					if q < '!' || q > '~' {
						t.Fatalf("read %d: quality byte %d outside Phred+33", i, q)
					}
				}
				if data.Headers[i] == "" {
					t.Fatalf("read %d has an empty header", i)
				}
			}
			if spec.Dist != "" && data.Stats.MinReadLength == data.Stats.MaxReadLength {
				t.Errorf("all %d reads have length %d; want a %s distribution", c.Reads, data.Stats.MinReadLength, spec.Dist)
			}
			if want := ComputeReadStats(data.Sequences, data.QualityScores); data.Stats != want {
				t.Errorf("Stats = %+v, want recomputation %+v", data.Stats, want)
			}
			if data.Stats.GCContent < 0.4 || data.Stats.GCContent > 0.6 {
				t.Errorf("GC content %.3f: the engine draws uniform bases, want about 0.5", data.Stats.GCContent)
			}
		})
	}
}

// TestGenerateSyntheticData_SeedDeterminism checks that the same seed gives
// identical data and a different seed different data.
func TestGenerateSyntheticData_SeedDeterminism(t *testing.T) {
	f := NewE2EFramework()
	for _, p := range modeledPlatforms {
		a := generate(t, f, p, 42, 50)
		b := generate(t, f, p, 42, 50)
		if !reflect.DeepEqual(a, b) {
			t.Errorf("%s: same seed gave different data", p)
		}
		c := generate(t, f, p, 43, 50)
		if reflect.DeepEqual(a.Sequences, c.Sequences) {
			t.Errorf("%s: seeds 42 and 43 gave identical sequences", p)
		}
	}
}

// fastqTextStats recomputes mean quality, GC content and base count
// directly from FASTQ text, independently of parseFASTQ.
func fastqTextStats(t *testing.T, path string) (meanQ, gc float64, bases int64, reads int) {
	t.Helper()
	fh, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()
	sc := bufio.NewScanner(fh)
	sc.Buffer(make([]byte, 1<<20), 1<<26)
	var qsum, qn, gcN, acgt int64
	for line := 0; sc.Scan(); line++ {
		text := sc.Text()
		switch line % 4 {
		case 1:
			reads++
			bases += int64(len(text))
			acgt += int64(strings.Count(text, "A") + strings.Count(text, "T"))
			n := int64(strings.Count(text, "G") + strings.Count(text, "C"))
			gcN += n
			acgt += n
		case 3:
			for i := 0; i < len(text); i++ {
				qsum += int64(text[i]) - 33
				qn++
			}
		}
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return float64(qsum) / float64(qn), float64(gcN) / float64(acgt), bases, reads
}

// TestStats_MatchFASTQText checks that the reported mean quality, GC content
// and base count equal a recomputation from the FASTQ text written for the
// same reads.
func TestStats_MatchFASTQText(t *testing.T) {
	f := NewE2EFramework()
	for _, p := range modeledPlatforms {
		data := generate(t, f, p, 7, 100)
		dir := t.TempDir()
		r1, r2 := filepath.Join(dir, "r1.fastq"), filepath.Join(dir, "r2.fastq")
		if err := f.WriteFASTQ(data, r1, r2); err != nil {
			t.Fatal(err)
		}
		meanQ, gc, bases, reads := fastqTextStats(t, r1)
		if reads != data.Stats.Reads || bases != data.Stats.Bases {
			t.Errorf("%s: FASTQ has %d reads / %d bases, Stats %d / %d", p, reads, bases, data.Stats.Reads, data.Stats.Bases)
		}
		if math.Abs(meanQ-data.Stats.MeanQuality) > 1e-9 {
			t.Errorf("%s: MeanQuality %.6f, FASTQ text gives %.6f", p, data.Stats.MeanQuality, meanQ)
		}
		if math.Abs(gc-data.Stats.GCContent) > 1e-9 {
			t.Errorf("%s: GCContent %.6f, FASTQ text gives %.6f", p, data.Stats.GCContent, gc)
		}
		b1, _ := os.ReadFile(r1)
		b2, _ := os.ReadFile(r2)
		if !bytes.Equal(b1, b2) {
			t.Errorf("%s: R2 differs from R1; WriteFASTQ documents them as identical", p)
		}

		loaded, err := LoadFASTQ(r1, 0)
		if err != nil {
			t.Fatal(err)
		}
		if !reflect.DeepEqual(loaded.Sequences, data.Sequences) || !reflect.DeepEqual(loaded.QualityScores, data.QualityScores) ||
			!reflect.DeepEqual(loaded.Headers, data.Headers) || loaded.Stats != data.Stats {
			t.Errorf("%s: LoadFASTQ(WriteFASTQ(data)) differs from data", p)
		}
		if part, err := LoadFASTQ(r1, 10); err != nil || len(part.Sequences) != 10 {
			t.Errorf("%s: LoadFASTQ(max 10) = %v reads, err %v", p, len(part.Sequences), err)
		}
	}
}

func TestComputeReadStats(t *testing.T) {
	seqs := [][]byte{[]byte("GGGG"), []byte("ATATAT"), []byte("ACGTNNNNNN")}
	quals := [][]byte{[]byte("IIII"), []byte("!!!!!!"), []byte("5555555555")} // Q40, Q0, Q20
	s := ComputeReadStats(seqs, quals)
	want := ReadStats{
		Reads: 3, Bases: 20, MeanReadLength: 20.0 / 3, MinReadLength: 4, MaxReadLength: 10,
		N50:         10, // lengths 10, 6, 4: the longest read already holds half of the 20 bases
		MeanQuality: (4*40 + 6*0 + 10*20) / 20.0,
		GCContent:   (4 + 0 + 2) / 14.0,
	}
	if s != want {
		t.Errorf("ComputeReadStats = %+v, want %+v", s, want)
	}
	if got := ComputeReadStats([][]byte{[]byte("AAA"), []byte("A"), []byte("AAAAAA")}, nil).N50; got != 6 {
		t.Errorf("N50 of lengths 3, 1, 6 (total 10) = %d, want 6", got)
	}
	if got := ComputeReadStats(nil, nil); got != (ReadStats{}) {
		t.Errorf("empty input: %+v, want zero", got)
	}
}

func TestParseFASTQ_Malformed(t *testing.T) {
	cases := map[string]string{
		"no @":       "read\nACGT\n+\nIIII\n",
		"truncated":  "@r\nACGT\n+\n",
		"length":     "@r\nACGT\n+\nIII\n",
		"no plus":    "@r\nACGT\n-\nIIII\n",
		"second bad": "@r\nA\n+\nI\nx\nA\n+\nI\n",
	}
	for name, text := range cases {
		if _, err := parseFASTQ(strings.NewReader(text), 0); err == nil {
			t.Errorf("%s: want error", name)
		}
	}
}

// TestRunBenchmark checks that the benchmark reports measured values that
// agree with each other and with the generated reads.
func TestRunBenchmark(t *testing.T) {
	f := NewE2EFramework()
	for _, p := range modeledPlatforms {
		c := f.CreatePlatformConfig(p, "CEU")
		b, err := f.RunBenchmark(c)
		if err != nil {
			t.Fatalf("%s: %v", p, err)
		}
		data, err := f.GenerateSyntheticData(c)
		if err != nil {
			t.Fatal(err)
		}
		if b.Stats != data.Stats {
			t.Errorf("%s: benchmark stats %+v differ from the same config's data %+v", p, b.Stats, data.Stats)
		}
		if b.WallTime <= 0 || b.AllocatedBytes == 0 {
			t.Errorf("%s: WallTime %v, AllocatedBytes %d; want both positive", p, b.WallTime, b.AllocatedBytes)
		}
		if want := float64(b.Stats.Bases) / b.WallTime.Seconds(); math.Abs(b.BasesPerSecond-want) > 1e-6*want {
			t.Errorf("%s: BasesPerSecond %.1f, want Bases/WallTime %.1f", p, b.BasesPerSecond, want)
		}
		// Regression: CPU use was once read from a metric of CPU time
		// available (GOMAXPROCS × wall), reporting ~GOMAXPROCS × 100%.
		if b.CPUSeconds < 0 || b.CPUPercent < 0 || b.CPUPercent > 100 {
			t.Errorf("%s: CPUSeconds %.3f, CPUPercent %.1f; want >= 0 and <= 100%%", p, b.CPUSeconds, b.CPUPercent)
		}
	}
	if _, err := f.RunBenchmark(nil); err == nil {
		t.Error("RunBenchmark(nil): want error")
	}
}

func TestCPUSecondsSinceStart_Monotonic(t *testing.T) {
	before := cpuSecondsSinceStart()
	deadline := time.Now().Add(50 * time.Millisecond)
	for time.Now().Before(deadline) {
	}
	if after := cpuSecondsSinceStart(); runtime.GOOS != "windows" && after <= before {
		t.Errorf("CPU time did not advance during a busy loop: before=%.4f after=%.4f", before, after)
	}
}

func TestMLClassifier_EmptyIsError(t *testing.T) {
	f := NewE2EFramework()
	data := generate(t, f, IlluminaNovaSeq6000, 3, 20)
	c := f.CreateMLClassifier()
	if err := c.Train(data, &SyntheticData{}); err == nil {
		t.Error("Train with an empty set: want error")
	}
	if _, err := c.TestAccuracy(&SyntheticData{}, data); err == nil {
		t.Error("TestAccuracy with an empty set: want error, not a made-up 0.5")
	}
	if _, _, err := HoldoutAccuracy(data, &SyntheticData{Sequences: [][]byte{[]byte("A")}}); err == nil {
		t.Error("HoldoutAccuracy with one read: want error")
	}
}

// allGC replaces every A with G and T with C.
func allGC(d *SyntheticData) *SyntheticData {
	out := *d
	out.Sequences = make([][]byte, len(d.Sequences))
	for i, s := range d.Sequences {
		out.Sequences[i] = bytes.Map(func(r rune) rune {
			switch r {
			case 'A':
				return 'G'
			case 'T':
				return 'C'
			}
			return r
		}, s)
	}
	return &out
}

// TestMLClassifier_TestAccuracy_UsesTrainedModel checks that a trained
// classifier scores held-out data with the fitted model.
func TestMLClassifier_TestAccuracy_UsesTrainedModel(t *testing.T) {
	f := NewE2EFramework()
	trainA, trainB := generate(t, f, IlluminaNovaSeq6000, 1, 100), generate(t, f, IlluminaNovaSeq6000, 2, 100)
	testA, testB := generate(t, f, IlluminaNovaSeq6000, 3, 100), generate(t, f, IlluminaNovaSeq6000, 4, 100)

	c := f.CreateMLClassifier()
	if err := c.Train(trainA, allGC(trainB)); err != nil {
		t.Fatal(err)
	}
	if acc, err := c.TestAccuracy(testA, allGC(testB)); err != nil || acc < 0.99 {
		t.Errorf("held-out accuracy = %.3f (err %v), want >= 0.99", acc, err)
	}
	// Swapping the labels must invert the result, which cross-validation on
	// the held-out set (ignoring the fitted model) would not do.
	if acc, err := c.TestAccuracy(allGC(testB), testA); err != nil || acc > 0.01 {
		t.Errorf("label-swapped held-out accuracy = %.3f (err %v), want <= 0.01", acc, err)
	}
	// Untrained: 5-fold cross-validation.
	if acc, err := f.CreateMLClassifier().TestAccuracy(testA, allGC(testB)); err != nil || acc < 0.99 {
		t.Errorf("cross-validated accuracy = %.3f (err %v), want >= 0.99", acc, err)
	}
}

// TestGCSkewed checks the positive-control transform.
func TestGCSkewed(t *testing.T) {
	data := generate(t, NewE2EFramework(), IlluminaNovaSeq6000, 5, 500)
	orig := make([][]byte, len(data.Sequences))
	for i, s := range data.Sequences {
		orig[i] = append([]byte(nil), s...)
	}

	for _, target := range []float64{0.6, 0.7, 0.9} {
		skewed := GCSkewed(data, target, 7)
		if math.Abs(skewed.Stats.GCContent-target) > 0.01 {
			t.Errorf("target %.2f: GC %.4f", target, skewed.Stats.GCContent)
		}
		if skewed.Stats != ComputeReadStats(skewed.Sequences, skewed.QualityScores) {
			t.Errorf("target %.2f: Stats not recomputed", target)
		}
	}
	for i := range orig {
		if !bytes.Equal(data.Sequences[i], orig[i]) {
			t.Fatalf("GCSkewed modified input sequence %d", i)
		}
	}
	if a, b := GCSkewed(data, 0.7, 7), GCSkewed(data, 0.7, 7); !reflect.DeepEqual(a.Sequences, b.Sequences) {
		t.Error("GCSkewed is not deterministic for a fixed seed")
	}
	if same := GCSkewed(data, data.Stats.GCContent-0.1, 7); !reflect.DeepEqual(same.Sequences, orig) {
		t.Error("a target below the current GC content changed the sequences")
	}
}

// TestGCSkewControl_Power is the power check for classifier comparisons. Two
// independent seeded draws of 300 bp MiSeq reads (1000 per set, 1000
// held-out reads) are compared after skewing the second set's GC content to
// p. At p = 0.6 (a 0.1 shift; the per-read GC standard deviation at 300 bp
// is 0.029) the classifier must exceed 0.9 accuracy. At p = 0.5 there is
// no skew, so accuracy must lie within 4σ of 0.5 (σ = sqrt(0.25/1000)).
func TestGCSkewControl_Power(t *testing.T) {
	f := NewE2EFramework()
	a := generate(t, f, IlluminaMiSeq, 101, 1000)
	b := generate(t, f, IlluminaMiSeq, 202, 1000)

	acc, n, err := HoldoutAccuracy(a, GCSkewed(b, 0.6, 9))
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("p=0.6: accuracy %.4f on %d held-out reads", acc, n)
	if acc <= 0.9 {
		t.Errorf("p=0.6: accuracy %.4f, want > 0.9", acc)
	}

	acc, n, err = HoldoutAccuracy(a, GCSkewed(b, 0.5, 9))
	if err != nil {
		t.Fatal(err)
	}
	band := ChanceBand(n, 4)
	t.Logf("p=0.5: accuracy %.4f on %d held-out reads (band 0.5 ± %.4f)", acc, n, band)
	if math.Abs(acc-0.5) > band {
		t.Errorf("p=0.5: accuracy %.4f outside 0.5 ± %.4f", acc, band)
	}
}

func TestChanceBand(t *testing.T) {
	if got, want := ChanceBand(2000, 4), 4*math.Sqrt(0.25/2000); math.Abs(got-want) > 1e-15 {
		t.Errorf("ChanceBand(2000, 4) = %v, want %v", got, want)
	}
	if !math.IsInf(ChanceBand(0, 4), 1) {
		t.Error("ChanceBand(0) should be +Inf")
	}
}

// TestWriteVCF checks that every site's REF is the reference base at POS.
func TestWriteVCF(t *testing.T) {
	f := NewE2EFramework()
	data := generate(t, f, IlluminaNovaSeq6000, 11, 0)
	dir := t.TempDir()
	ref, vcf := filepath.Join(dir, "ref.fa"), filepath.Join(dir, "v.vcf")
	if err := f.WriteReferenceFromData(data, ref); err != nil {
		t.Fatal(err)
	}
	if err := f.WriteVCF(data, vcf); err != nil {
		t.Fatal(err)
	}
	refText, _ := os.ReadFile(ref)
	contig := strings.Split(string(refText), "\n")[1]
	vcfText, _ := os.ReadFile(vcf)

	sites := 0
	for _, line := range strings.Split(strings.TrimSpace(string(vcfText)), "\n") {
		if strings.HasPrefix(line, "#") {
			if strings.HasPrefix(line, "##contig=") && line != "##contig=<ID=chr1,length="+strconv.Itoa(len(contig))+">" {
				t.Errorf("contig line %q, want length %d", line, len(contig))
			}
			continue
		}
		fields := strings.Split(line, "\t")
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			t.Fatal(err)
		}
		if got := string(contig[pos-1]); fields[3] != got {
			t.Errorf("POS %d: REF %s, reference base %s", pos, fields[3], got)
		}
		if alt, _ := transition(fields[3][0]); fields[4] != string(alt) {
			t.Errorf("POS %d: ALT %s is not the transition of %s", pos, fields[4], fields[3])
		}
		sites++
	}
	if want := min(vcfMaxSites, len(contig)/vcfSiteSpacing); sites != want {
		t.Errorf("%d sites, want %d", sites, want)
	}
}

func TestWriteBAM(t *testing.T) {
	f := NewE2EFramework()
	data := generate(t, f, IlluminaNovaSeq6000, DefaultSeed, 0)
	fname := t.TempDir() + "/test.bam"
	err := f.WriteBAM(data, fname)
	if _, lookErr := exec.LookPath("samtools"); lookErr != nil {
		if !errors.Is(err, ErrBAMUnsupported) {
			t.Fatalf("WriteBAM without samtools: err = %v, want ErrBAMUnsupported", err)
		}
		t.Skip("samtools not installed; BAM conversion not exercised")
	}
	if err != nil {
		t.Fatalf("WriteBAM error: %v", err)
	}

	// The BAM must hold one record per read; the old header stub held none.
	out, err := exec.Command("samtools", "view", "-c", fname).Output()
	if err != nil {
		t.Fatalf("samtools view -c: %v", err)
	}
	if got := strings.TrimSpace(string(out)); got != strconv.Itoa(len(data.Sequences)) {
		t.Errorf("BAM record count = %s, want %d", got, len(data.Sequences))
	}
	if err := exec.Command("samtools", "quickcheck", fname).Run(); err != nil {
		t.Errorf("samtools quickcheck failed: %v", err)
	}
}

// TestWriteBAM_NoSamtools checks that WriteBAM refuses rather than writing a
// fake file when samtools is missing.
func TestWriteBAM_NoSamtools(t *testing.T) {
	orig := lookPath
	t.Cleanup(func() { lookPath = orig })
	lookPath = func(string) (string, error) { return "", exec.ErrNotFound }

	f := NewE2EFramework()
	data := generate(t, f, IlluminaNovaSeq6000, DefaultSeed, 10)
	fname := t.TempDir() + "/test.bam"
	if err := f.WriteBAM(data, fname); !errors.Is(err, ErrBAMUnsupported) {
		t.Fatalf("WriteBAM err = %v, want ErrBAMUnsupported", err)
	}
	if _, err := os.Stat(fname); !os.IsNotExist(err) {
		t.Errorf("WriteBAM wrote %s despite missing samtools", fname)
	}
}

// TestWriteSAM checks the header and that every read appears as an unmapped
// record with its sequence and qualities.
func TestWriteSAM(t *testing.T) {
	f := NewE2EFramework()
	for _, p := range []SequencingPlatform{IlluminaNovaSeq6000, ONTPromethION24} {
		data := generate(t, f, p, DefaultSeed, 50)
		fname := t.TempDir() + "/test.sam"
		if err := f.WriteSAM(data, fname); err != nil {
			t.Fatalf("WriteSAM error: %v", err)
		}
		content, err := os.ReadFile(fname)
		if err != nil {
			t.Fatal(err)
		}

		var header, records []string
		for _, line := range strings.Split(strings.TrimRight(string(content), "\n"), "\n") {
			if strings.HasPrefix(line, "@") {
				header = append(header, line)
			} else {
				records = append(records, line)
			}
		}
		if len(header) < 3 || !strings.HasPrefix(header[0], "@HD\tVN:") {
			t.Fatalf("missing header lines: %q", header)
		}
		if want := "@SQ\tSN:chr1\tLN:" + strconv.FormatInt(data.Stats.Bases, 10); header[1] != want {
			t.Errorf("@SQ line = %q, want %q", header[1], want)
		}
		if want := "@RG\tID:" + samReadGroupID + "\tSM:CEU\tPL:" + samPlatform(p); header[2] != want {
			t.Errorf("@RG line = %q, want %q", header[2], want)
		}
		if len(records) != len(data.Sequences) {
			t.Fatalf("records = %d, want one per read (%d)", len(records), len(data.Sequences))
		}
		for i, rec := range records {
			fields := strings.Split(rec, "\t")
			if len(fields) != 12 {
				t.Fatalf("record %d has %d fields, want 11 mandatory plus RG: %q", i, len(fields), rec)
			}
			if fields[11] != "RG:Z:"+samReadGroupID {
				t.Errorf("record %d tag = %q, want RG:Z:%s", i, fields[11], samReadGroupID)
			}
			if want := strings.Fields(data.Headers[i])[0]; fields[0] != want {
				t.Errorf("record %d QNAME = %q, want %q", i, fields[0], want)
			}
			if got := strings.Join(fields[1:9], " "); got != "4 * 0 0 * * 0 0" {
				t.Errorf("record %d fields 2-9 = %q, want unmapped \"4 * 0 0 * * 0 0\"", i, got)
			}
			if fields[9] != string(data.Sequences[i]) || fields[10] != string(data.QualityScores[i]) {
				t.Errorf("record %d SEQ/QUAL do not match read %d", i, i)
			}
		}
	}
	if err := f.WriteSAM(&SyntheticData{}, t.TempDir()+"/x.sam"); err == nil {
		t.Error("WriteSAM with no reads should fail")
	}
}
