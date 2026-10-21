//go:build !release

package testing

import (
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

// phred encodes Phred values as Phred+33.
func phred(q ...int) []byte {
	b := make([]byte, len(q))
	for i, v := range q {
		b[i] = byte(v + 33)
	}
	return b
}

// feature returns the named element of QualityFeatures(seq, qual).
func feature(t *testing.T, f []float64, name string) float64 {
	t.Helper()
	for i, n := range qualityFeatureNames {
		if n == name {
			return f[i]
		}
	}
	t.Fatalf("no quality feature %q", name)
	return 0
}

// TestQualityFeatures_KnownRead checks every kind of quality feature on a
// read worked out by hand.
func TestQualityFeatures_KnownRead(t *testing.T) {
	seq := []byte("ACGTN")
	qual := phred(37, 37, 11, 11, 2)
	f := QualityFeatures(seq, qual)
	if len(f) != len(QualityFeatureNames()) {
		t.Fatalf("%d features, %d names", len(f), len(QualityFeatureNames()))
	}
	mean := (37 + 37 + 11 + 11 + 2) / 5.0
	sd := math.Sqrt((37*37+37*37+11*11+11*11+2*2)/5.0 - mean*mean)
	want := map[string]float64{
		"qual_frac_Q37": 0.4, "qual_frac_Q11": 0.4, "qual_frac_Q2": 0.2, "qual_frac_Q25": 0,
		// 5 cycles in 10 blocks: cycle i is in block 2i; odd blocks are empty
		// and take the read mean.
		"mean_qual_block_1/10": 37, "mean_qual_block_2/10": mean, "mean_qual_block_3/10": 37,
		"mean_qual_block_5/10": 11, "mean_qual_block_7/10": 11, "mean_qual_block_9/10": 2,
		"mean_qual": mean, "sd_qual": sd,
		"lowq_runs_per_100bp(Q<20)": 20, "lowq_longest_run(Q<20)": 3, "lowq_tail_run(Q<20)": 3,
		"qual_change_frac": 0.5, "n_frac": 0.2,
		"mean_qual_A": 37, "mean_qual_C": 37, "mean_qual_G": 11, "mean_qual_T": 11,
	}
	for name, w := range want {
		if got := feature(t, f, name); math.Abs(got-w) > 1e-12 {
			t.Errorf("%s = %v, want %v", name, got, w)
		}
	}

	// Two low-quality runs, the read ends on high quality, no N, a base
	// that does not occur takes the mean, and Q>=50 shares the last bin.
	f = QualityFeatures([]byte("AAAAAA"), phred(5, 30, 5, 5, 30, 60))
	for name, w := range map[string]float64{
		"lowq_runs_per_100bp(Q<20)": 100 * 2 / 6.0, "lowq_longest_run(Q<20)": 2, "lowq_tail_run(Q<20)": 0,
		"n_frac": 0, "mean_qual_C": (5 + 30 + 5 + 5 + 30 + 60) / 6.0, "qual_frac_Q>=50": 1 / 6.0,
	} {
		if got := feature(t, f, name); math.Abs(got-w) > 1e-12 {
			t.Errorf("second read: %s = %v, want %v", name, got, w)
		}
	}
	if f := QualityFeatures(nil, nil); len(f) != len(qualityFeatureNames) {
		t.Errorf("empty read: %d features", len(f))
	}
}

// TestFeatureNames checks that every set's names match its vectors and
// that 4-mer names follow kmerIndex.
func TestFeatureNames(t *testing.T) {
	seq, qual := []byte("ACGTACGGTTAC"), phred(37, 37, 37, 25, 25, 11, 37, 37, 37, 2, 37, 37)
	for _, set := range AllFeatureSets {
		names, err := FeatureNames(set)
		if err != nil {
			t.Fatal(err)
		}
		x, err := ReadFeatures(set, seq, qual)
		if err != nil {
			t.Fatal(err)
		}
		if len(names) != len(x) {
			t.Errorf("%s: %d names for %d features", set, len(names), len(x))
		}
	}
	for _, k := range []string{"AAAA", "ACGT", "TTTT", "GATC"} {
		if got := kmerNames[kmerIndex([]byte(k))]; got != k {
			t.Errorf("kmerNames[kmerIndex(%s)] = %s", k, got)
		}
	}
	if _, err := FeatureNames("bogus"); err == nil {
		t.Error("unknown feature set: want error")
	}
	if _, err := ReadFeatures(FeatureQuality, seq, nil); err == nil {
		t.Error("quality features without qualities: want error")
	}
	if _, err := ReadFeatures(FeatureSequence, seq, nil); err != nil {
		t.Errorf("sequence features without qualities: %v", err)
	}
}

// gaussianSets returns n rows per class of d features, N(0,1) except that
// class 1's first feature has mean shift.
func gaussianSets(n, d int, shift float64, seed uint64) ([][]float64, []int) {
	rng := rand.New(rand.NewPCG(seed, 1))
	var x [][]float64
	var y []int
	for c := range 2 {
		for range n {
			row := make([]float64, d)
			for k := range row {
				row[k] = rng.NormFloat64()
			}
			if c == 1 {
				row[0] += shift
			}
			row[d-1] = 7 // constant: must be dropped, not divide by zero
			x, y = append(x, row), append(y, c)
		}
	}
	return x, y
}

func accuracy(c BinaryClassifier, x [][]float64, y []int) float64 {
	ok := 0
	for i := range x {
		if c.Predict(x[i]) == y[i] {
			ok++
		}
	}
	return float64(ok) / float64(len(x))
}

// TestLogisticRegression checks power on separable data (a 4-SD shift:
// Bayes accuracy Φ(2) = 0.977), chance on identical distributions,
// determinism and input validation.
func TestLogisticRegression(t *testing.T) {
	xTr, yTr := gaussianSets(1000, 20, 4, 1)
	xTe, yTe := gaussianSets(1000, 20, 4, 2)
	lr := NewLogisticRegression()
	if err := lr.Fit(xTr, yTr); err != nil {
		t.Fatal(err)
	}
	if acc := accuracy(lr, xTe, yTe); acc < 0.96 {
		t.Errorf("separable: held-out accuracy %.4f, want >= 0.96 (Bayes 0.977)", acc)
	}
	w := lr.Weights()
	for k := 1; k < len(w); k++ {
		if math.Abs(w[k]) >= math.Abs(w[0]) || math.IsNaN(w[k]) {
			t.Errorf("weight %d = %v not below the informative weight %v", k, w[k], w[0])
		}
	}
	if w[len(w)-1] != 0 {
		t.Errorf("constant feature weight %v, want 0", w[len(w)-1])
	}
	again := NewLogisticRegression()
	_ = again.Fit(xTr, yTr)
	if !reflect.DeepEqual(again.Weights(), w) {
		t.Error("two fits on the same data differ")
	}

	xTr, yTr = gaussianSets(1000, 20, 0, 3)
	xTe, yTe = gaussianSets(1000, 20, 0, 4)
	lr = NewLogisticRegression()
	if err := lr.Fit(xTr, yTr); err != nil {
		t.Fatal(err)
	}
	if acc := accuracy(lr, xTe, yTe); math.Abs(acc-0.5) > ChanceBand(len(xTe), 4) {
		t.Errorf("inseparable: accuracy %.4f outside 0.5 ± %.4f", acc, ChanceBand(len(xTe), 4))
	}

	if p := NewLogisticRegression().PredictProba([]float64{1}); p != 0.5 {
		t.Errorf("untrained PredictProba = %v", p)
	}
	for _, bad := range []struct {
		x [][]float64
		y []int
	}{
		{nil, nil},
		{[][]float64{{1}}, []int{0, 1}},
		{[][]float64{{1}, {2}}, []int{0, 0}},
		{[][]float64{{1}, {2}}, []int{0, 2}},
		{[][]float64{{1}, {2, 3}}, []int{0, 1}},
		{[][]float64{{}, {}}, []int{0, 1}},
	} {
		if err := NewLogisticRegression().Fit(bad.x, bad.y); err == nil {
			t.Errorf("Fit(%v, %v): want error", bad.x, bad.y)
		}
	}
}

// stackedRows builds combined vectors: 4-mer counts of a random read with
// GC probability gc, then three continuous features with mean qmean.
func stackedRows(rng *rand.Rand, n int, gc, qmean float64, label int) ([][]float64, []int) {
	var x [][]float64
	var y []int
	for range n {
		var sb strings.Builder
		for range 100 {
			if rng.Float64() < gc {
				sb.WriteByte("GC"[rng.IntN(2)])
			} else {
				sb.WriteByte("AT"[rng.IntN(2)])
			}
		}
		row := Extract4merFrequencies(sb.String())
		for range 3 {
			row = append(row, qmean+rng.NormFloat64())
		}
		x, y = append(x, row), append(y, label)
	}
	return x, y
}

// TestStackedClassifier checks the combined classifier sees a difference
// in either part, and nothing when there is none.
func TestStackedClassifier(t *testing.T) {
	cases := []struct {
		name         string
		gc, qmean    float64
		separable    bool
		minAccuracy  float64
		inChanceBand bool
	}{
		{"sequence differs", 0.7, 0, true, 0.95, false},
		{"quality differs", 0.5, 3, true, 0.95, false},
		{"neither differs", 0.5, 0, false, 0, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			rng := rand.New(rand.NewPCG(5, 6))
			build := func() ([][]float64, []int) {
				x0, y0 := stackedRows(rng, 600, 0.5, 0, 0)
				x1, y1 := stackedRows(rng, 600, c.gc, c.qmean, 1)
				return append(x0, x1...), append(y0, y1...)
			}
			xTr, yTr := build()
			xTe, yTe := build()
			s := NewStackedClassifier()
			if err := s.Fit(xTr, yTr); err != nil {
				t.Fatal(err)
			}
			acc := accuracy(s, xTe, yTe)
			if c.separable && acc < c.minAccuracy {
				t.Errorf("accuracy %.4f, want >= %.2f", acc, c.minAccuracy)
			}
			if c.inChanceBand && math.Abs(acc-0.5) > ChanceBand(len(xTe), 4) {
				t.Errorf("accuracy %.4f outside 0.5 ± %.4f", acc, ChanceBand(len(xTe), 4))
			}
		})
	}
	if err := NewStackedClassifier().Fit([][]float64{make([]float64, kmerVocabSize)}, []int{0}); err == nil {
		t.Error("combined vector without continuous features: want error")
	}
}

// TestHoldoutAccuracyFor_Calibration is the calibration of every feature
// set on generated NovaSeq reads: two independent draws are within the
// chance band, as are two halves of one draw (the real-versus-real null
// applied to generated reads); a GC skew is seen by the sets with
// sequence features, the other quality model by the sets with quality
// features, and a one-Phred quality shift by the quality sets.
func TestHoldoutAccuracyFor_Calibration(t *testing.T) {
	f := NewE2EFramework()
	ca := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	ca.Reads, ca.Seed = 1000, 11
	cb := ca.clone()
	cb.Seed = 12
	a, err := f.GenerateSyntheticData(ca)
	if err != nil {
		t.Fatal(err)
	}
	b, err := f.GenerateSyntheticData(cb)
	if err != nil {
		t.Fatal(err)
	}
	alt, err := f.GenerateSyntheticData(cb.AlternativeQualityModel())
	if err != nil {
		t.Fatal(err)
	}
	even, odd := SplitAlternate(a)
	if len(even.QualityScores) != 500 || len(odd.QualityScores) != 500 || even.Stats.Reads != 500 {
		t.Fatalf("SplitAlternate: %d/%d qualities, stats %+v", len(even.QualityScores), len(odd.QualityScores), even.Stats)
	}
	gc := GCSkewed(b, 0.7, 3)
	shifted := QualityShifted(b, 1)

	for _, set := range AllFeatureSets {
		null, err := HoldoutAccuracyFor(set, a, b, 3)
		if err != nil {
			t.Fatal(err)
		}
		if band := ChanceBand(null.Holdout, 4); null.Holdout != 1000 || math.Abs(null.Accuracy-0.5) > band {
			t.Errorf("%s null: accuracy %.4f on %d held out, want 1000 within 0.5 ± %.4f", set, null.Accuracy, null.Holdout, band)
		}
		if len(null.Top) != 3 || null.Classifier != ClassifierName(set) {
			t.Errorf("%s null: %d top features, classifier %q", set, len(null.Top), null.Classifier)
		}
		halves, err := HoldoutAccuracyFor(set, even, odd, 0)
		if err != nil {
			t.Fatal(err)
		}
		if band := ChanceBand(halves.Holdout, 4); math.Abs(halves.Accuracy-0.5) > band {
			t.Errorf("%s halves of one draw: accuracy %.4f outside 0.5 ± %.4f", set, halves.Accuracy, band)
		}

		gcAcc, err := HoldoutAccuracyFor(set, a, gc, 1)
		if err != nil {
			t.Fatal(err)
		}
		qAcc, err := HoldoutAccuracyFor(set, a, alt, 1)
		if err != nil {
			t.Fatal(err)
		}
		shAcc, err := HoldoutAccuracyFor(set, a, shifted, 0)
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("%s: null %.4f, halves %.4f, GC %.4f, quality model %.4f, Q-1 %.4f", set, null.Accuracy, halves.Accuracy, gcAcc.Accuracy, qAcc.Accuracy, shAcc.Accuracy)
		if set != FeatureQuality && gcAcc.Accuracy <= 0.9 {
			t.Errorf("%s: GC control %.4f, want > 0.9", set, gcAcc.Accuracy)
		}
		if set != FeatureSequence {
			if qAcc.Accuracy <= 0.9 || shAcc.Accuracy <= 0.9 {
				t.Errorf("%s: quality controls %.4f / %.4f, want > 0.9", set, qAcc.Accuracy, shAcc.Accuracy)
			}
		} else if math.Abs(shAcc.Accuracy-0.5) > ChanceBand(shAcc.Holdout, 4) {
			t.Errorf("sequence set sees a quality-only shift: %.4f", shAcc.Accuracy)
		}
		// The top sequence feature of the GC control is a GC-rich 4-mer
		// that is more frequent in the skewed set.
		if set == FeatureSequence {
			top := gcAcc.Top[0]
			if strings.Count(top.Name, "G")+strings.Count(top.Name, "C") < 3 || top.EffectSize <= 0 || top.MeanB <= top.MeanA {
				t.Errorf("GC control top feature %+v, want a GC-rich 4-mer enriched in the skewed set", top)
			}
		}
		if set == FeatureQuality && !strings.HasPrefix(qAcc.Top[0].Name, "qual_frac_") && !strings.HasPrefix(qAcc.Top[0].Name, "qual_change") {
			t.Errorf("quality-model control top feature %+v, want a quality-value feature", qAcc.Top[0])
		}
	}

	noQual := &SyntheticData{Sequences: a.Sequences}
	if _, err := HoldoutAccuracyFor(FeatureQuality, noQual, b, 0); err == nil {
		t.Error("quality set on reads without qualities: want error")
	}
}

// TestAlternativeQualityModel checks the quality control configs.
func TestAlternativeQualityModel(t *testing.T) {
	f := NewE2EFramework()
	il := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	if alt := il.AlternativeQualityModel(); alt.IlluminaQualityModel != "linear" || il.IlluminaQualityModel != "" {
		t.Errorf("illumina: %q (original %q)", alt.IlluminaQualityModel, il.IlluminaQualityModel)
	}
	if back := il.AlternativeQualityModel().AlternativeQualityModel(); back.IlluminaQualityModel != "binned" {
		t.Errorf("linear -> %q, want binned", back.IlluminaQualityModel)
	}
	ont := f.CreatePlatformConfig(ONTPromethION24, "CEU")
	if alt := ont.AlternativeQualityModel(); alt.ErrorTier != "" || ont.ErrorTier != "sup" {
		t.Errorf("ont: %q (original %q)", alt.ErrorTier, ont.ErrorTier)
	}
	if f.CreatePlatformConfig(PacBioRevio, "CEU").AlternativeQualityModel() != nil {
		t.Error("unmodeled platform has an alternative quality model")
	}
}

// revcomp returns the reverse complement of s.
func revcomp(s string) string {
	b := make([]byte, len(s))
	for i := range s {
		b[len(s)-1-i] = map[byte]byte{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[s[i]]
	}
	return string(b)
}

// readContig returns the contig of fasta whose forward or reverse strand
// holds read with at most maxMismatch mismatches (N never matches), or ""
// when none does. Brute force: the fixture is small.
func readContig(read string, contigs map[string]string, maxMismatch int) string {
	for name, seq := range contigs {
		for _, s := range []string{seq, revcomp(seq)} {
			for p := 0; p+len(read) <= len(s); p++ {
				mm := 0
				for i := 0; i < len(read) && mm <= maxMismatch; i++ {
					if read[i] != s[p+i] {
						mm++
					}
				}
				if mm <= maxMismatch {
					return name
				}
			}
		}
	}
	return ""
}

// loadFixture writes a two-contig reference and returns its path and
// sequences.
func loadFixture(t *testing.T) (string, map[string]string) {
	t.Helper()
	path := filepath.Join(t.TempDir(), "ref.fa")
	if err := WriteRandomReference(path, []int{3000, 2000}, 42); err != nil {
		t.Fatal(err)
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	contigs := map[string]string{}
	var name string
	for _, line := range strings.Split(string(raw), "\n") {
		if strings.HasPrefix(line, ">") {
			name = line[1:]
			continue
		}
		contigs[name] += line
	}
	return path, contigs
}

// TestGenerateSyntheticData_Reference checks generation from a reference:
// every read comes from the fixture (within the error model's
// substitutions), a region restricts reads to it, and uniform reads do not
// match the fixture.
func TestGenerateSyntheticData_Reference(t *testing.T) {
	path, contigs := loadFixture(t)
	if len(contigs) != 2 || len(contigs["ctg1"]) != 3000 || len(contigs["ctg2"]) != 2000 {
		t.Fatalf("fixture contigs %d", len(contigs))
	}
	f := NewE2EFramework()
	c := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	c.Reads, c.Seed = 60, 7
	c.ReferencePath = path
	if c.BaseSource() != BasesReference {
		t.Errorf("BaseSource = %s", c.BaseSource())
	}
	data, err := f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatal(err)
	}
	seen := map[string]int{}
	for i, s := range data.Sequences {
		name := readContig(string(s), contigs, 15)
		if name == "" {
			t.Fatalf("read %d does not come from the fixture: %s", i, s)
		}
		seen[name]++
	}
	if seen["ctg1"] == 0 || seen["ctg2"] == 0 {
		t.Errorf("reads per contig %v, want both contigs sampled", seen)
	}

	c.ReferenceRegion = "ctg2:101-1900"
	data, err = f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatal(err)
	}
	region := map[string]string{"ctg2": contigs["ctg2"][100:1900]}
	for i, s := range data.Sequences {
		if readContig(string(s), region, 15) == "" {
			t.Fatalf("region read %d is not from ctg2:101-1900", i)
		}
	}

	c.ReferencePath, c.ReferenceRegion = "", ""
	if c.BaseSource() != BasesUniform {
		t.Errorf("BaseSource = %s", c.BaseSource())
	}
	uniform, err := f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatal(err)
	}
	if readContig(string(uniform.Sequences[0]), contigs, 15) != "" {
		t.Error("a uniform random read matches the fixture")
	}

	c.ReferencePath = filepath.Join(t.TempDir(), "missing.fa")
	if _, err := f.GenerateSyntheticData(c); err == nil {
		t.Error("missing reference: want error")
	}
}
