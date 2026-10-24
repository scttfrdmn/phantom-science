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

// pairFeature returns the named element of a PairFeatures vector.
func pairFeature(t *testing.T, f []float64, name string) float64 {
	t.Helper()
	for i, n := range pairFeatureNames {
		if n == name {
			return f[i]
		}
	}
	t.Fatalf("no pair feature %q", name)
	return 0
}

func TestPairFeatures_KnownPair(t *testing.T) {
	seq := []byte("ACGTACGTAC")
	// Read 1: mean 30, 2 of 10 bases below Q20; read 2: mean 24, 5 below.
	q1 := phred(37, 37, 37, 11, 11, 37, 37, 37, 25, 31)
	q2 := phred(37, 11, 11, 11, 11, 11, 37, 37, 37, 37)
	f, err := PairFeatures(seq, q1, seq, q2)
	if err != nil {
		t.Fatal(err)
	}
	if len(f) != len(PairFeatureNames()) || len(f) != 2*len(qualityFeatureNames)+len(pairOnlyFeatureNames) {
		t.Fatalf("%d features, want %d", len(f), len(PairFeatureNames()))
	}
	if !reflect.DeepEqual(f[:len(qualityFeatureNames)], QualityFeatures(seq, q1)) ||
		!reflect.DeepEqual(f[len(qualityFeatureNames):2*len(qualityFeatureNames)], QualityFeatures(seq, q2)) {
		t.Error("the first features are not the two mates' quality features")
	}
	for name, want := range map[string]float64{
		"r1:mean_qual":                 30,
		"r2:mean_qual":                 24,
		"r1:qual_frac_Q11":             0.2,
		"r2:qual_frac_Q11":             0.5,
		"r1:qual_change_frac":          4.0 / 9,
		"r2:qual_change_frac":          2.0 / 9,
		"r2:lowq_longest_run(Q<20)":    5,
		"r1:mean_qual_block_9/10":      25,
		"r2:mean_qual_block_2/10":      11,
		"pair_mean_qual_diff(R2-R1)":   -6,
		"pair_abs_mean_qual_diff":      6,
		"pair_lowq_frac_product(Q<20)": 0.2 * 0.5,
		"pair_lowq_frac_min(Q<20)":     0.2,
	} {
		if got := pairFeature(t, f, name); math.Abs(got-want) > 1e-9 {
			t.Errorf("%s = %g, want %g", name, got, want)
		}
	}
	// Mates whose block means move together correlate positively.
	if c := pairFeature(t, f, "pair_block_mean_corr"); c <= 0 || c > 1 {
		t.Errorf("block correlation %g, want in (0, 1]", c)
	}
	same, _ := PairFeatures(seq, q1, seq, q1)
	if c := pairFeature(t, same, "pair_block_mean_corr"); math.Abs(c-1) > 1e-12 {
		t.Errorf("identical mates: block correlation %g, want 1", c)
	}
	flat, _ := PairFeatures(seq, phred(37, 37, 37, 37, 37, 37, 37, 37, 37, 37), seq, q2)
	if c := pairFeature(t, flat, "pair_block_mean_corr"); c != 0 {
		t.Errorf("constant mate: block correlation %g, want 0", c)
	}

	if _, err := PairFeatures(seq, q1[:5], seq, q2); err == nil {
		t.Error("PairFeatures accepted a mate without qualities")
	}
	if _, err := ReadFeatures(FeaturePair, seq, q1); err == nil {
		t.Error("ReadFeatures accepted the pair set for one read")
	}
	names, err := FeatureNames(FeaturePair)
	if err != nil || !reflect.DeepEqual(names, PairFeatureNames()) {
		t.Errorf("FeatureNames(pair) = %v, %v", names, err)
	}
}

func TestPairName(t *testing.T) {
	for h, want := range map[string]string{
		"DRR276626.1 A00718:197:HGT2NDSXY:1:1101:1118:1078/1": "DRR276626.1",
		"DRR394916.1 1/2":                "DRR394916.1",
		"read_7/1":                       "read_7",
		"read_7/2\tRG:Z:x":               "read_7",
		"A00718:1:FC:1:1:1:1 1:N:0:ACGT": "A00718:1:FC:1:1:1:1",
		"x/3":                            "x/3",
	} {
		if got := pairName(h); got != want {
			t.Errorf("pairName(%q) = %q, want %q", h, got, want)
		}
	}
}

// writeFASTQText writes FASTQ records with the given names to a file.
func writeFASTQText(t *testing.T, path string, names ...string) {
	t.Helper()
	var b strings.Builder
	for _, n := range names {
		b.WriteString("@" + n + "\nACGTACGTAC\n+\nFFFF,FFFF:\n")
	}
	if err := os.WriteFile(path, []byte(b.String()), 0600); err != nil {
		t.Fatal(err)
	}
}

func TestLoadPairedFASTQ(t *testing.T) {
	dir := t.TempDir()
	p1, p2 := filepath.Join(dir, "r1.fastq"), filepath.Join(dir, "r2.fastq")
	writeFASTQText(t, p1, "a/1", "b/1", "c/1")
	writeFASTQText(t, p2, "a/2", "b/2", "c/2")
	d, err := LoadPairedFASTQ(p1, p2, 0)
	if err != nil {
		t.Fatal(err)
	}
	if !d.Paired() || len(d.Sequences) != 3 || d.Mates.Headers[2] != "c/2" || d.Mates.Stats.Reads != 3 {
		t.Fatalf("pairs %+v", d)
	}
	if d, err := LoadPairedFASTQ(p1, p2, 2); err != nil || len(d.Mates.Sequences) != 2 {
		t.Errorf("maxPairs 2: %v, %v", d, err)
	}

	short := filepath.Join(dir, "short.fastq")
	writeFASTQText(t, short, "a/2", "b/2")
	if _, err := LoadPairedFASTQ(p1, short, 0); err == nil || !strings.Contains(err.Error(), "3 reads 1 but 2 reads 2") {
		t.Errorf("unequal files: %v", err)
	}
	other := filepath.Join(dir, "other.fastq")
	writeFASTQText(t, other, "a/2", "x/2", "c/2")
	if _, err := LoadPairedFASTQ(p1, other, 0); err == nil || !strings.Contains(err.Error(), "record 2") {
		t.Errorf("mismatched names: %v", err)
	}
	if _, err := LoadPairedFASTQ(p1, filepath.Join(dir, "missing.fastq"), 0); err == nil {
		t.Error("missing read-2 file accepted")
	}
}

// TestGenerateSyntheticData_Paired generates read pairs from a reference
// fixture and checks that they are pairs: one mate per read, matching
// names, read 1 forward and read 2 reverse-complemented from the same
// fragment, reproducible from the seed, and written to two files.
func TestGenerateSyntheticData_Paired(t *testing.T) {
	dir := t.TempDir()
	ref := filepath.Join(dir, "ref.fa")
	if err := WriteRandomReference(ref, []int{50000}, 5); err != nil {
		t.Fatal(err)
	}
	genome := readFASTASequence(t, ref)
	f := NewE2EFramework()
	c := f.CreatePlatformConfig(IlluminaNovaSeq6000, "CEU")
	c.Seed, c.Reads, c.ReadLength, c.PairedEnd, c.ReferencePath = 11, 300, 151, true, ref
	d, err := f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatal(err)
	}
	if !d.Paired() || len(d.Sequences) != 300 || d.Mates.Stats.Reads != 300 || d.Stats.Reads != 300 {
		t.Fatalf("paired data: %d reads, mates %v", len(d.Sequences), d.Mates != nil)
	}
	if d.Mates.Platform != IlluminaNovaSeq6000 || d.Mates.Seed != 11 {
		t.Errorf("mates labelled %s seed %d", d.Mates.Platform, d.Mates.Seed)
	}
	t.Logf("first pair headers %q / %q", d.Headers[0], d.Mates.Headers[0])
	// With a 151 bp read and few errors, the first 20 bases of most reads
	// occur in the genome, and the mates of a pair on opposite strands:
	// read 1 forward and read 2 reverse-complemented, or the other way round
	// for fragments from the reverse strand. Each mate is checked on its
	// first 20 cycles: the last cycles have the most errors, and in the
	// default profile (DRR394916) about a quarter of reads have an error in
	// their last 20.
	in := func(b []byte) bool { return strings.Contains(genome, string(b[:20])) }
	inRC := func(b []byte) bool { return strings.Contains(genome, reverseComplement(string(b[:20]))) }
	opposite := 0
	for i := range d.Sequences {
		r1, r2 := d.Sequences[i], d.Mates.Sequences[i]
		if (in(r1) && inRC(r2)) || (inRC(r1) && in(r2)) {
			opposite++
		}
	}
	if opposite < 250 {
		t.Errorf("of 300 pairs, %d have mates on opposite strands of the reference; want most", opposite)
	}

	again, err := f.GenerateSyntheticData(c)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(again.Mates.Sequences, d.Mates.Sequences) || !reflect.DeepEqual(again.QualityScores, d.QualityScores) {
		t.Error("paired generation is not reproducible from the seed")
	}

	p1, p2 := filepath.Join(dir, "g1.fastq"), filepath.Join(dir, "g2.fastq")
	if err := f.WriteFASTQ(d, p1, p2); err != nil {
		t.Fatal(err)
	}
	back, err := LoadPairedFASTQ(p1, p2, 0)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(back.Mates.Sequences, d.Mates.Sequences) || !reflect.DeepEqual(back.Sequences, d.Sequences) {
		t.Error("written pairs do not read back")
	}
}

// readFASTASequence returns the concatenated sequence lines of a FASTA.
func readFASTASequence(t *testing.T, path string) string {
	t.Helper()
	b, err := os.ReadFile(path) // #nosec G304 -- test fixture
	if err != nil {
		t.Fatal(err)
	}
	var s strings.Builder
	for _, line := range strings.Split(string(b), "\n") {
		if !strings.HasPrefix(line, ">") {
			s.WriteString(strings.ToUpper(strings.TrimSpace(line)))
		}
	}
	return s.String()
}

func reverseComplement(s string) string {
	b := make([]byte, len(s))
	for i := range s {
		c := s[len(s)-1-i]
		switch c {
		case 'A':
			c = 'T'
		case 'C':
			c = 'G'
		case 'G':
			c = 'C'
		case 'T':
			c = 'A'
		}
		b[i] = c
	}
	return string(b)
}

// coupledPairs returns n read pairs whose mates share a quality class:
// both good (few Q11 calls) or both poor (many), each with probability 1/2.
// With shuffle the read-2 set is permuted, so each mate's distribution is
// unchanged but the classes of the mates are independent.
func coupledPairs(n int, seed uint64, shuffle bool) *SyntheticData {
	rng := rand.New(rand.NewPCG(seed, 2)) // #nosec G404 -- deterministic test fixture
	read := func(poor bool) ([]byte, []byte) {
		s, q := make([]byte, 100), make([]byte, 100)
		p := 0.01
		if poor {
			p = 0.15
		}
		for i := range s {
			s[i] = "ACGT"[rng.IntN(4)]
			q[i] = 'F'
			if rng.Float64() < p {
				q[i] = ','
			}
		}
		return s, q
	}
	r1, r2 := &SyntheticData{}, &SyntheticData{}
	for range n {
		poor := rng.IntN(2) == 0
		s, q := read(poor)
		r1.Sequences, r1.QualityScores = append(r1.Sequences, s), append(r1.QualityScores, q)
		s, q = read(poor)
		r2.Sequences, r2.QualityScores = append(r2.Sequences, s), append(r2.QualityScores, q)
	}
	if shuffle {
		rng.Shuffle(n, func(i, j int) {
			r2.Sequences[i], r2.Sequences[j] = r2.Sequences[j], r2.Sequences[i]
			r2.QualityScores[i], r2.QualityScores[j] = r2.QualityScores[j], r2.QualityScores[i]
		})
	}
	r1.Mates = r2
	return r1
}

// TestHoldoutAccuracyFor_Pair checks the pair feature set: it is at chance
// on two draws of the same pairs, it sees that the mates' qualities go
// together (coupled against shuffled mates, which the read-1 quality set
// cannot see), and it needs paired sets.
func TestHoldoutAccuracyFor_Pair(t *testing.T) {
	const n = 4000
	band := ChanceBand(n, 4)
	null, err := HoldoutAccuracyFor(FeaturePair, coupledPairs(n, 1, false), coupledPairs(n, 2, false), 3)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(null.Accuracy-0.5) > band || null.Holdout != n || null.Classifier != "logistic regression" || len(null.Top) != 3 {
		t.Errorf("null %+v, want 0.5 ± %.4f on %d held out", null, band, n)
	}
	coupled, shuffled := coupledPairs(n, 3, false), coupledPairs(n, 4, true)
	r1, err := HoldoutAccuracyFor(FeatureQuality, coupled, shuffled, 0)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(r1.Accuracy-0.5) > band {
		t.Errorf("read-1 quality set on coupled vs shuffled mates: %.4f, want chance (read 1 is the same)", r1.Accuracy)
	}
	p, err := HoldoutAccuracyFor(FeaturePair, coupled, shuffled, 5)
	if err != nil {
		t.Fatal(err)
	}
	if p.Accuracy < 0.6 {
		t.Errorf("pair set on coupled vs shuffled mates: %.4f, want > 0.6", p.Accuracy)
	}
	if !strings.HasPrefix(p.Top[0].Name, "pair_") {
		t.Errorf("top separating feature %q, want a pair_ feature", p.Top[0].Name)
	}
	single := &SyntheticData{Sequences: coupled.Sequences, QualityScores: coupled.QualityScores}
	if _, err := HoldoutAccuracyFor(FeaturePair, coupled, single, 0); err == nil {
		t.Error("pair set accepted an unpaired set")
	}

	// The controls carry the mates.
	even, odd := SplitAlternate(coupled)
	if !even.Paired() || !odd.Paired() || even.Mates.Stats.Reads != n/2 {
		t.Error("SplitAlternate dropped the mates")
	}
	if q := QualityShifted(coupled, 1); !q.Paired() || q.Mates.QualityScores[0][0] != coupled.Mates.QualityScores[0][0]-1 {
		t.Error("QualityShifted did not shift the mates")
	}
	if g := GCSkewed(coupled, 0.8, 1); !g.Paired() || g.Mates.Stats.GCContent < 0.75 {
		t.Errorf("GCSkewed did not skew the mates (GC %.3f)", g.Mates.Stats.GCContent)
	}
}

func TestTrimmedTo(t *testing.T) {
	d := coupledPairs(10, 1, false)
	d.Stats = ComputeReadStats(d.Sequences, d.QualityScores)
	tr := TrimmedTo(d, 60)
	if tr.Stats.MaxReadLength != 60 || tr.Mates.Stats.MaxReadLength != 60 || len(tr.QualityScores[3]) != 60 {
		t.Errorf("trimmed stats %+v / %+v", tr.Stats, tr.Mates.Stats)
	}
	if string(tr.Sequences[0]) != string(d.Sequences[0][:60]) || d.Stats.MaxReadLength != 100 {
		t.Error("TrimmedTo changed the input or kept the wrong bases")
	}
	if same := TrimmedTo(d, 200); same.Stats != d.Stats {
		t.Errorf("trimming to more than the read length changed the reads: %+v", same.Stats)
	}
}
