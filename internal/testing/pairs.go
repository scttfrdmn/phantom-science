package testing

import (
	"fmt"
	"math"
	"strings"
)

// Paired reads. A paired read set is a SyntheticData of the first reads
// with Mates holding the second reads in the same order: Mates.Sequences[i]
// is the mate of Sequences[i]. GenerateSyntheticData produces one when
// PlatformConfig.PairedEnd is set, LoadPairedFASTQ from two FASTQ files.

// pairOnlyFeatureNames name the features of PairFeatures after the two
// mates' quality features.
var pairOnlyFeatureNames = []string{
	"pair_mean_qual_diff(R2-R1)",
	"pair_abs_mean_qual_diff",
	fmt.Sprintf("pair_lowq_frac_product(Q<%d)", lowQualThreshold),
	fmt.Sprintf("pair_lowq_frac_min(Q<%d)", lowQualThreshold),
	"pair_block_mean_corr",
}

// pairFeatureNames is the name of each element of PairFeatures.
var pairFeatureNames = func() []string {
	var n []string
	for _, m := range []string{"r1:", "r2:"} {
		for _, q := range qualityFeatureNames {
			n = append(n, m+q)
		}
	}
	return append(n, pairOnlyFeatureNames...)
}()

// PairFeatureNames returns the names of the PairFeatures elements.
func PairFeatureNames() []string { return append([]string(nil), pairFeatureNames...) }

// PairFeatures returns the quality profile of a read pair: QualityFeatures
// of read 1, QualityFeatures of read 2, then
//
//   - read 2's mean quality minus read 1's, and its absolute value;
//   - the product and the minimum of the mates' low-quality (Phred < 20)
//     base fractions, which are large only when both mates are poor;
//   - the Pearson correlation of the mates' 10 cycle-block mean qualities
//     (0 when either mate's block means are constant).
//
// A logistic regression on the two mates' features alone is additive in
// them; the last five features let it see how the mates' qualities go
// together (a poor cluster gives two poor reads). None of them depends on
// the bases or on the insert size.
func PairFeatures(seq1, qual1, seq2, qual2 []byte) ([]float64, error) {
	if len(qual1) != len(seq1) || len(qual2) != len(seq2) {
		return nil, fmt.Errorf("feature set %s needs qualities: mates have %d/%d and %d/%d bases/qualities",
			FeaturePair, len(seq1), len(qual1), len(seq2), len(qual2))
	}
	f1, f2 := QualityFeatures(seq1, qual1), QualityFeatures(seq2, qual2)
	f := make([]float64, 0, len(pairFeatureNames))
	f = append(append(f, f1...), f2...)
	meanIdx := qualHistBins + qualBlocks
	diff := f2[meanIdx] - f1[meanIdx]
	low1, low2 := lowQualFraction(qual1), lowQualFraction(qual2)
	f = append(f, diff, math.Abs(diff), low1*low2, math.Min(low1, low2),
		pearson(f1[qualHistBins:qualHistBins+qualBlocks], f2[qualHistBins:qualHistBins+qualBlocks]))
	return f, nil
}

// lowQualFraction is the fraction of qualities below lowQualThreshold.
func lowQualFraction(qual []byte) float64 {
	if len(qual) == 0 {
		return 0
	}
	n := 0
	for _, c := range qual {
		if int(c)-33 < lowQualThreshold {
			n++
		}
	}
	return float64(n) / float64(len(qual))
}

// pearson returns the correlation of x and y, or 0 when either is constant.
func pearson(x, y []float64) float64 {
	n := float64(len(x))
	if n == 0 || len(y) != len(x) {
		return 0
	}
	var mx, my float64
	for i := range x {
		mx += x[i]
		my += y[i]
	}
	mx /= n
	my /= n
	var sxy, sxx, syy float64
	for i := range x {
		dx, dy := x[i]-mx, y[i]-my
		sxy += dx * dy
		sxx += dx * dx
		syy += dy * dy
	}
	if sxx < 1e-12 || syy < 1e-12 {
		return 0
	}
	return sxy / math.Sqrt(sxx*syy)
}

// pairName returns the read name shared by both mates of a FASTQ header
// (without '@'): the first whitespace-separated token without a trailing
// "/1" or "/2".
func pairName(header string) string {
	name, _, _ := strings.Cut(header, " ")
	name, _, _ = strings.Cut(name, "\t")
	if n := len(name); n > 2 && name[n-2] == '/' && (name[n-1] == '1' || name[n-1] == '2') {
		name = name[:n-2]
	}
	return name
}

// LoadPairedFASTQ reads up to maxPairs read pairs (all when maxPairs <= 0)
// from the FASTQ files of read 1 and read 2 and returns the read-1 set
// with the read-2 set as its Mates. The files must hold the same number of
// records (within maxPairs) with matching read names.
func LoadPairedFASTQ(path1, path2 string, maxPairs int) (*SyntheticData, error) {
	r1, err := LoadFASTQ(path1, maxPairs)
	if err != nil {
		return nil, fmt.Errorf("read 1: %w", err)
	}
	r2, err := LoadFASTQ(path2, maxPairs)
	if err != nil {
		return nil, fmt.Errorf("read 2: %w", err)
	}
	if err := pairUp(r1, r2); err != nil {
		return nil, fmt.Errorf("%s and %s: %w", path1, path2, err)
	}
	return r1, nil
}

// pairUp checks that r2 holds the mates of r1's reads, by count and name,
// and sets r1.Mates.
func pairUp(r1, r2 *SyntheticData) error {
	if len(r1.Sequences) != len(r2.Sequences) {
		return fmt.Errorf("%d reads 1 but %d reads 2", len(r1.Sequences), len(r2.Sequences))
	}
	for i := range r1.Headers {
		if i < len(r2.Headers) && pairName(r1.Headers[i]) != pairName(r2.Headers[i]) {
			return fmt.Errorf("record %d: read names %q and %q differ", i+1, pairName(r1.Headers[i]), pairName(r2.Headers[i]))
		}
	}
	r1.Mates = r2
	return nil
}

// deinterleave splits interleaved paired FASTQ records (read 1, read 2,
// read 1, ...) into a read-1 set with the read-2 set as its Mates.
func deinterleave(d *SyntheticData) (*SyntheticData, error) {
	if len(d.Sequences)%2 != 0 {
		return nil, fmt.Errorf("interleaved FASTQ holds an odd number of records (%d)", len(d.Sequences))
	}
	r1, r2 := split(d)
	r1.Stats = ComputeReadStats(r1.Sequences, r1.QualityScores)
	r2.Stats = ComputeReadStats(r2.Sequences, r2.QualityScores)
	if err := pairUp(r1, r2); err != nil {
		return nil, err
	}
	return r1, nil
}

// Paired reports whether d holds read pairs (Mates set, one mate per read).
func (d *SyntheticData) Paired() bool {
	return d != nil && d.Mates != nil && len(d.Mates.Sequences) == len(d.Sequences)
}

// TrimmedTo returns a copy of d (and of its Mates) with every read longer
// than length cut to its first length bases and qualities; Stats are
// recomputed. Comparing two runs whose read lengths differ (2×150 and
// 2×151) on their common cycles keeps the length itself, and the last
// cycle's qualities, out of the comparison.
func TrimmedTo(d *SyntheticData, length int) *SyntheticData {
	out := *d
	out.Sequences = make([][]byte, len(d.Sequences))
	out.QualityScores = make([][]byte, len(d.QualityScores))
	for i, s := range d.Sequences {
		out.Sequences[i] = s[:min(len(s), length)]
	}
	for i, q := range d.QualityScores {
		out.QualityScores[i] = q[:min(len(q), length)]
	}
	out.Stats = ComputeReadStats(out.Sequences, out.QualityScores)
	if d.Mates != nil {
		out.Mates = TrimmedTo(d.Mates, length)
	}
	return &out
}
