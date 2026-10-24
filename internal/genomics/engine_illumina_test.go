//go:build !release

package genomics

// Tests of the Illumina read model (#27): binned qualities, quality-matched
// substitution errors, insert sizes, and paired reads that are the two ends
// of one fragment.

import (
	"bytes"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/illumina"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// binnedConfig is goldenConfig with the default Illumina model (binned
// qualities, quality errors) and paired-end reads.
func binnedConfig(dir string) types.SequencerConfig {
	c := goldenConfig(dir)
	c.Illumina = types.IlluminaSpec{}
	c.PairedEnd = true
	c.ReadLength = 150
	c.Seed = 27
	return c
}

// binChars are the NovaSeq 6000 quality values (Q2, Q11, Q25, Q37) in Phred+33.
const binChars = "#,:F"

// mateFiles generates c to files and returns the R1 and R2 records.
func mateFiles(t *testing.T, c types.SequencerConfig) (r1, r2 []fastqRec) {
	t.Helper()
	for name, data := range fileOutputs(t, c) {
		switch {
		case strings.Contains(name, "_R1_"):
			r1 = append(r1, parseFASTQ(t, string(data))...)
		case strings.Contains(name, "_R2_"):
			r2 = append(r2, parseFASTQ(t, string(data))...)
		}
	}
	if len(r1) == 0 || len(r1) != len(r2) {
		t.Fatalf("%d R1 and %d R2 reads", len(r1), len(r2))
	}
	return r1, r2
}

// The default Illumina output uses only the four bins, writes every Q2 base
// as N and no other base as N, degrades with cycle, and R2 is worse than R1.
// An HMM model's qualities are replaced by the binned model too.
func TestIllumina_DefaultBinnedQualities(t *testing.T) {
	for name, edit := range map[string]func(*types.SequencerConfig){
		"uniform": func(*types.SequencerConfig) {},
		"hmm":     func(c *types.SequencerConfig) { c.HMMModelPath = trainedModel(t) },
	} {
		t.Run(name, func(t *testing.T) {
			c := binnedConfig(t.TempDir())
			c.ReadsPerBatch, c.TotalBatches = 2000, 1
			edit(&c)
			r1, r2 := mateFiles(t, c)
			meanQ := func(recs []fastqRec, lo, hi int) float64 {
				var s, n float64
				for _, r := range recs {
					for _, q := range r.qual[lo:hi] {
						s += float64(q - '!')
						n++
					}
				}
				return s / n
			}
			for _, recs := range [][]fastqRec{r1, r2} {
				for _, r := range recs {
					for i := range r.qual {
						if strings.IndexByte(binChars, r.qual[i]) < 0 {
							t.Fatalf("quality %q is not a NovaSeq 6000 bin", r.qual[i])
						}
						if (r.qual[i] == '#') != (r.seq[i] == 'N') {
							t.Fatalf("base %c with quality %c: Q2 must be exactly the N calls", r.seq[i], r.qual[i])
						}
					}
				}
				if early, late := meanQ(recs, 0, 30), meanQ(recs, 120, 150); !(early > late+0.5) {
					t.Errorf("mean quality cycles 1-30 %.2f, 121-150 %.2f: want a decrease", early, late)
				}
			}
			// DRR276626 (NovaSeq 6000) measured R1 36.22 and R2 36.08.
			m1, m2 := meanQ(r1, 0, 150), meanQ(r2, 0, 150)
			if math.Abs(m1-36.22) > 0.15 || math.Abs(m2-36.08) > 0.15 || !(m2 < m1) {
				t.Errorf("mean quality R1 %.2f, R2 %.2f; want 36.22 and 36.08 (±0.15), R2 below R1", m1, m2)
			}
			t.Logf("mean quality R1 %.2f, R2 %.2f", m1, m2)
		})
	}
}

// The empirical substitution count per quality bin matches the model's
// fitted, cycle-dependent rates (illumina.SubstitutionRate × BaseFactor of
// the true base), within 99 % of a Poisson-binomial. Truth is the error-free
// run with the same seed: paired file output draws each fragment from the
// pair stream, so turning errors off changes only the mate streams
// (qualities and errors), not the fragments.
func TestIllumina_ErrorRateMatchesQuality(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.ReadsPerBatch, c.TotalBatches = 20000, 1
	o1, o2 := mateFiles(t, c)
	c.Illumina.Errors = types.IlluminaErrorsNone
	t1, t2 := mateFiles(t, c)
	baseIdx := map[byte]int{'A': 0, 'C': 1, 'G': 2, 'T': 3}
	var bases, subs [256]int
	var mean, variance [256]float64
	for mate, pair := range [][2][]fastqRec{{o1, t1}, {o2, t2}} {
		obs, truth := pair[0], pair[1]
		for i := range obs {
			if strings.Contains(truth[i].seq, "N") {
				t.Fatalf("mate %d: error-free read has N", mate+1)
			}
			for j := range obs[i].seq {
				qc := obs[i].qual[j]
				bases[qc]++
				if obs[i].seq[j] != truth[i].seq[j] {
					subs[qc]++
				}
				if q := qc - '!'; q != illumina.QNoCall {
					p := illumina.SubstitutionRate(q, j) * illumina.NovaSeqErrors.BaseFactor[baseIdx[truth[i].seq[j]]]
					mean[qc] += p
					variance[qc] += p * (1 - p)
				}
			}
		}
	}
	var total, totalSubs, noCalls int
	for _, b := range []byte(binChars) {
		q := b - '!'
		n, k := bases[b], subs[b]
		total += n
		if q == illumina.QNoCall {
			noCalls = n
			if k != n {
				t.Errorf("Q2: %d of %d bases changed, want all (N)", k, n)
			}
			continue
		}
		totalSubs += k
		ci := 2.5758 * math.Sqrt(variance[b])
		if math.Abs(float64(k)-mean[b]) > ci {
			t.Errorf("Q%d: %d substitutions over %d bases, want %.1f ± %.1f", q, k, n, mean[b], ci)
		}
		t.Logf("Q%d: %d bases, %d substitutions (expected %.1f ± %.1f)", q, n, k, mean[b], ci)
	}
	t.Logf("overall: substitutions %.4f%%, N calls %.5f%% of %d bases", 100*float64(totalSubs)/float64(total), 100*float64(noCalls)/float64(total), total)
}

// recordingSampler is a FragmentSampler over uniform bases that records the
// fragments it returns.
type recordingSampler struct{ frags [][]byte }

func (s *recordingSampler) sample(r *rand.Rand, n int) (string, int64, bool, []byte, error) {
	f := randomBases(r, n)
	s.frags = append(s.frags, f)
	return "chrT", int64(len(s.frags)) * 1000, false, bytes.Clone(f), nil
}

// Each mate's binned qualities are drawn for its own bases in sequencing
// orientation (#35): R1 for the start of the fragment, R2 for the reverse
// complement of its end, from the mate's stream, before errors. Errors off,
// the pipeline's mate gives exactly illumina.ParamsFor(mate).QualitiesForClass
// on those bases with the same stream and class.
func TestIllumina_QualitiesFollowMateBases(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.Illumina.Errors = types.IlluminaErrorsNone
	e := NewEngine(c)
	defer e.Close()
	p, err := e.newPipeline(GenerateOpts{}, nil)
	if err != nil {
		t.Fatal(err)
	}
	src := rand.New(rand.NewPCG(35, 1)) // #nosec G404 -- test data
	const acgt = "ACGT"
	for i := range 50 {
		frag := make([]byte, 200+i)
		for j := range frag {
			frag[j] = acgt[src.IntN(4)]
		}
		for _, mate := range []int{1, 2} {
			want := append([]byte(nil), frag[:c.ReadLength]...)
			if mate == 2 {
				want = seqread.ReverseComplement(frag[len(frag)-c.ReadLength:])
			}
			seed := uint64(100*i + mate)
			class := i%7 - 1                                                                                 // -1: drawn from the mate's stream
			seq, qual := p.mate(rand.New(rand.NewPCG(seed, 2)), frag, c.ReadLength, mate, class)             // #nosec G404 -- test stream
			wantQ := illumina.ParamsFor(mate).QualitiesForClass(rand.New(rand.NewPCG(seed, 2)), want, class) // #nosec G404 -- test stream
			if string(seq) != string(want) || string(qual) != string(wantQ) {
				t.Fatalf("fragment %d mate %d: qualities not drawn for the mate's bases in sequencing orientation", i, mate)
			}
		}
	}
}

// With errors off, R1 is the start of the fragment and R2 the reverse
// complement of its end; insert sizes follow the configured distribution.
// The fragment hook (FragmentSampler) shows the truth.
func TestIllumina_MatesAreFragmentEnds(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 100, 5000, 1
	c.Illumina.Errors = types.IlluminaErrorsNone
	c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 300, SD: 40}
	e := NewEngine(c)
	defer e.Close()
	s := &recordingSampler{}
	e.fragments = s.sample
	var buf bytes.Buffer
	if err := e.GenerateToWriter(&buf, GenerateOpts{}); err != nil {
		t.Fatal(err)
	}
	recs := parseFASTQ(t, buf.String())
	if len(recs) != 2*len(s.frags) {
		t.Fatalf("%d reads for %d fragments", len(recs), len(s.frags))
	}
	L := c.ReadLength
	var sum, sq float64
	for i, f := range s.frags {
		r1, r2 := recs[2*i], recs[2*i+1]
		if r1.seq != string(f[:L]) {
			t.Fatalf("pair %d: R1 is not frag[:L]", i)
		}
		if got := string(seqread.ReverseComplement([]byte(r2.seq))); got != string(f[len(f)-L:]) {
			t.Fatalf("pair %d: revcomp(R2) is not frag[len-L:]", i)
		}
		if n1, n2 := strings.Fields(r1.header)[0], strings.Fields(r2.header)[0]; n1 != n2 {
			t.Fatalf("pair %d: mate names %s and %s differ", i, n1, n2)
		}
		if len(f) < L {
			t.Fatalf("fragment of %d bases shorter than the read", len(f))
		}
		sum += float64(len(f))
		sq += float64(len(f)) * float64(len(f))
	}
	n := float64(len(s.frags))
	m := sum / n
	sd := math.Sqrt(sq/n - m*m)
	if math.Abs(m-300) > 2.5758*40/math.Sqrt(n) || math.Abs(sd-40) > 2*2.5758*40/math.Sqrt(2*n) {
		t.Errorf("insert mean %.2f sd %.2f, want 300 ± %.2f and 40", m, sd, 2.5758*40/math.Sqrt(n))
	}
	t.Logf("insert sizes: mean %.2f sd %.2f over %d fragments", m, sd, len(s.frags))
}

// On the default (uniform-base) path, R1 and revcomp(R2) of a short
// fragment overlap consistently, in both file and stream output: some
// insert F in [L, 2L) makes R1[F−L:] equal the start of revcomp(R2).
func TestIllumina_MatesOverlapConsistently(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 250, 500, 1
	c.Illumina.Errors = types.IlluminaErrorsNone
	c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 350, SD: 50} // mostly overlapping mates
	check := func(t *testing.T, r1, r2 []fastqRec) {
		L := c.ReadLength
		overlapped := 0
		for i := range r1 {
			rc := string(seqread.ReverseComplement([]byte(r2[i].seq)))
			for f := L; f < 2*L; f++ {
				if r1[i].seq[f-L:] == rc[:2*L-f] {
					overlapped++
					break
				}
			}
		}
		// Inserts are 350 ± 50 with L = 250: P(F < 500) ≈ 1.
		if overlapped < len(r1)*99/100 {
			t.Errorf("%d of %d pairs overlap consistently", overlapped, len(r1))
		}
	}
	t.Run("files", func(t *testing.T) {
		r1, r2 := mateFiles(t, c)
		check(t, r1, r2)
	})
	t.Run("stream", func(t *testing.T) {
		recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{})))
		var r1, r2 []fastqRec
		for i := 0; i < len(recs); i += 2 {
			r1, r2 = append(r1, recs[i]), append(r2, recs[i+1])
		}
		check(t, r1, r2)
	})
}

// A positioned variant is applied once per fragment: when both mates cover
// it they carry the same allele. (Before #27 each mate drew its allele
// independently.)
func TestIllumina_MatesShareVariantAllele(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.ReadLength, c.ReadsPerBatch, c.TotalBatches = 250, 2000, 1
	c.Illumina.Errors = types.IlluminaErrorsNone
	c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 350, SD: 50} // fragments fit the test region
	// 16 bases, so that neither occurs by chance in the random bases.
	const ref, alt = "ACGTACGTACGTACGT", "TTTTTTTTTTTTTTTT"
	opts := GenerateOpts{
		RefRegion: &GenomicRegion{Chromosome: "chr1", Start: 1000, End: 1800},
		Variants:  []VariantRecord{{Chromosome: "chr1", Position: 1400, Ref: ref, Alt: alt, Freq: 0.5}},
	}
	recs := parseFASTQ(t, string(streamOutput(t, c, opts)))
	rcOf := func(s string) string { return string(seqread.ReverseComplement([]byte(s))) }
	allele := func(seq string) string {
		hasRef, hasAlt := strings.Contains(seq, ref), strings.Contains(seq, alt)
		switch {
		case hasRef && !hasAlt:
			return "ref"
		case hasAlt && !hasRef:
			return "alt"
		}
		return ""
	}
	counts := map[string]int{}
	for i := 0; i < len(recs); i += 2 {
		a1, a2 := allele(recs[i].seq), allele(rcOf(recs[i+1].seq))
		if a1 == "" || a2 == "" {
			continue
		}
		if a1 != a2 {
			t.Fatalf("pair %d: R1 carries %s, R2 %s", i/2, a1, a2)
		}
		counts[a1]++
	}
	if counts["ref"] < 20 || counts["alt"] < 20 {
		t.Errorf("pairs covering the variant with both mates: %v; want both alleles", counts)
	}
	t.Logf("pairs with both mates covering the variant: %v", counts)
}

// GenerateFASTQFile (the TUI path) writes exactly the file Generate writes
// for the same batch and mate, for both mates of binned paired output.
func TestIllumina_GenerateFASTQFileMatchesGenerate(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.TotalBatches = 2
	files := fileOutputs(t, c)
	for mate := 1; mate <= 2; mate++ {
		e := NewEngine(c)
		path := filepath.Join(c.OutputDir, "tui.fastq")
		if err := e.GenerateFASTQFile(path, 2, mate); err != nil {
			t.Fatal(err)
		}
		e.Close()
		got, err := os.ReadFile(path)
		if err != nil {
			t.Fatal(err)
		}
		name := filepath.Base(e.outputPath(2, mate))
		if want, ok := files[name]; !ok || !bytes.Equal(got, want) {
			t.Errorf("GenerateFASTQFile(batch 2, mate %d) differs from Generate's %s", mate, name)
		}
	}
}

// illuminaCommentHeader is the Illumina header line with the optional read
// group and read-ID fingerprint fields (#27).
var illuminaCommentHeader = regexp.MustCompile(`^@\S+ [12]:N:0:\d+(\tRG:Z:\S+)?(\tXC:Z:\S+)?(\tXD:Z:\S+)?$`)

// Read-ID fingerprint tags are SAM-style comment fields: the read name and
// the Illumina comment are unchanged, and the tags follow RG:Z:.
func TestIllumina_FingerprintTagsInComment(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.DatasetID = "ds 7"
	c.Enterprise.CustomerID = "CUST_001"
	c.Fingerprinting = types.FingerprintConfig{Enable: true, CustomerFingerprint: true, DatasetFingerprint: true, SteganographyLayers: []string{"read_ids"}}
	plain := binnedConfig(t.TempDir())
	for _, tc := range []struct {
		name   string
		recs   []fastqRec
		suffix string
	}{
		{"stream with RG", parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{SampleID: "s1"}))), "\tRG:Z:s1\tXC:Z:CUST_001\tXD:Z:ds_7"},
		{"files", func() []fastqRec { r1, r2 := mateFiles(t, c); return append(r1, r2...) }(), "\tXC:Z:CUST_001\tXD:Z:ds_7"},
		{"untagged", parseFASTQ(t, string(streamOutput(t, plain, GenerateOpts{}))), ""},
	} {
		for _, r := range tc.recs {
			if !illuminaCommentHeader.MatchString(r.header) || !strings.HasSuffix(r.header, tc.suffix) || strings.Contains(r.header, ":cust:") {
				t.Fatalf("%s: header %q, want the Illumina comment then %q", tc.name, r.header, tc.suffix)
			}
			if m := illuminaHeader.FindStringSubmatch(strings.Split(r.header, "\t")[0]); m == nil {
				t.Fatalf("%s: read name and comment %q changed", tc.name, r.header)
			}
		}
	}
}

// The illumina section is validated and only accepted on platform illumina.
func TestIllumina_ConfigChecks(t *testing.T) {
	for name, edit := range map[string]func(*types.SequencerConfig){
		"on ont":     func(c *types.SequencerConfig) { c.Platform = "ont"; c.PairedEnd = false },
		"bad mean":   func(c *types.SequencerConfig) { c.Illumina.InsertSize.Mean = -5 },
		"max < min":  func(c *types.SequencerConfig) { c.Illumina.InsertSize = types.InsertSizeSpec{Min: 400, Max: 300} },
		"bad errors": func(c *types.SequencerConfig) { c.Illumina.Errors = "indel" },
	} {
		c := goldenConfig(t.TempDir())
		edit(&c)
		e := NewEngine(c)
		if err := e.GenerateToWriter(&bytes.Buffer{}, GenerateOpts{}); err == nil {
			t.Errorf("%s: accepted", name)
		}
		e.Close()
	}
	// A region narrower than the fragments is an error, not a hang.
	c := binnedConfig(t.TempDir())
	e := NewEngine(c)
	defer e.Close()
	err := e.GenerateToWriter(&bytes.Buffer{}, GenerateOpts{RefRegion: &GenomicRegion{Chromosome: "chr1", Start: 0, End: 200}})
	if err == nil || !strings.Contains(err.Error(), "too small") {
		t.Errorf("region smaller than the insert: err = %v", err)
	}
}

// The two reads of a pair draw their quality classes together (#35): in
// file output (the classes come from the batch's pair stream, which the R1
// and R2 files each regenerate) and in interleaved streams, the mates'
// counts of bases below Q30 are correlated as in the fitted run (0.37;
// independent mates: 0).
func TestIllumina_MateQualitiesCoupled(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.ReadsPerBatch, c.TotalBatches = 3000, 2
	low := func(q string) float64 {
		n := 0.0
		for i := range len(q) {
			if q[i]-'!' < 30 {
				n++
			}
		}
		return n
	}
	check := func(t *testing.T, r1, r2 []fastqRec) {
		var sx, sy, sxx, syy, sxy float64
		for i := range r1 {
			x, y := low(r1[i].qual), low(r2[i].qual)
			sx, sy, sxx, syy, sxy = sx+x, sy+y, sxx+x*x, syy+y*y, sxy+x*y
		}
		n := float64(len(r1))
		got := (sxy - sx*sy/n) / math.Sqrt((sxx-sx*sx/n)*(syy-sy*sy/n))
		t.Logf("correlation of the mates' low-quality counts %.3f over %d pairs", got, len(r1))
		if got < 0.25 || got > 0.5 {
			t.Errorf("correlation of the mates' low-quality counts %.3f, want about 0.37", got)
		}
	}
	t.Run("files", func(t *testing.T) {
		r1, r2 := mateFiles(t, c)
		check(t, r1, r2)
	})
	t.Run("stream", func(t *testing.T) {
		recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{})))
		var r1, r2 []fastqRec
		for i := 0; i < len(recs); i += 2 {
			r1, r2 = append(r1, recs[i]), append(r2, recs[i+1])
		}
		check(t, r1, r2)
	})
}
