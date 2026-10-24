//go:build !release

package illumina

import (
	"math"
	"testing"
)

// The joint class weights are a distribution whose margins are the class
// weights of each mate, so coupling the mates leaves each mate's class
// mixture (and so its per-cycle fractions and mean quality) unchanged.
func TestPairs_MarginsAreClassWeights(t *testing.T) {
	for _, pr := range profiles {
		t.Run(pr.Name, func(t *testing.T) { testPairMargins(t, pr) })
	}
}

func testPairMargins(t *testing.T, pr *Profile) {
	pp, Read1, Read2 := pr.Pairs, pr.Read1, pr.Read2
	if !pp.Coupled() {
		t.Fatal("the mates' classes are not coupled")
	}
	if len(pp.Joint) != len(Read1.Classes) {
		t.Fatalf("%d rows, want %d (read 1 classes)", len(pp.Joint), len(Read1.Classes))
	}
	total := 0.0
	col := make([]float64, len(Read2.Classes))
	for i, row := range pp.Joint {
		if len(row) != len(Read2.Classes) {
			t.Fatalf("row %d: %d columns, want %d (read 2 classes)", i, len(row), len(Read2.Classes))
		}
		s := 0.0
		for j, w := range row {
			if !(w >= 0) {
				t.Errorf("Joint[%d][%d] = %v", i, j, w)
			}
			s += w
			col[j] += w
		}
		total += s
		if math.Abs(s-Read1.Classes[i].Weight) > 1e-5 {
			t.Errorf("read 1 class %d: margin %.6f, weight %.6f", i, s, Read1.Classes[i].Weight)
		}
	}
	for j, s := range col {
		if math.Abs(s-Read2.Classes[j].Weight) > 1e-5 {
			t.Errorf("read 2 class %d: margin %.6f, weight %.6f", j, s, Read2.Classes[j].Weight)
		}
	}
	if math.Abs(total-1) > 1e-5 {
		t.Errorf("joint weights sum to %.6f", total)
	}
}

// Classes maps uniforms onto the cells in proportion to their weights.
func TestPairs_ClassesFollowJoint(t *testing.T) {
	pp := PairParams{Joint: [][]float64{{0.5, 0.1}, {0.15, 0.25}}}
	for _, c := range []struct {
		u      float64
		k1, k2 int
	}{{0, 0, 0}, {0.49, 0, 0}, {0.5, 0, 1}, {0.59, 0, 1}, {0.61, 1, 0}, {0.76, 1, 1}, {0.999999, 1, 1}, {1, 1, 1}} {
		if k1, k2 := pp.Classes(c.u); k1 != c.k1 || k2 != c.k2 {
			t.Errorf("Classes(%v) = %d, %d; want %d, %d", c.u, k1, k2, c.k1, c.k2)
		}
	}
	if (PairParams{}).Coupled() {
		t.Error("the zero PairParams is coupled")
	}
}

// Pairs drawn with the joint classes have the measured coupling of the
// profile's run (#35): the correlation of the mates' Q2 + Q11 + Q25 counts
// and of their Q2 + Q11 counts (Phred < 20) within 0.03 of the measured
// (pair_summary.tsv), where independent mates have none; and each mate
// keeps its mean quality (within 4 standard errors of the independent
// draw).
func TestCalibration_MateCouplingMatchesNovaSeq(t *testing.T) {
	forProfiles(t, testMateCoupling)
}

func testMateCoupling(t *testing.T, pr *Profile) {
	want := map[string]float64{}
	for _, r := range readTSV(t, "pair_summary.tsv") {
		want[r["key"]] = num(t, r["value"])
	}
	loadReadCells(t)
	const pairs = 60000
	n := pr.Run.ReadLength
	type set struct {
		low30, low20 [2][]float64
		sum, sumSq   [2]float64 // of the reads' mean Q
	}
	draw := func(coupled bool) set {
		var s set
		r, sr := newRand(41), newRand(42)
		for range pairs {
			k := [2]int{-1, -1}
			if coupled {
				k[0], k[1] = pr.Pairs.Classes(r.Float64())
			}
			for m, p := range []QualityParams{pr.Read1, pr.Read2} {
				q := p.QualitiesForClass(r, humanBasesMate(sr, n, m), k[m])
				var a, b, sq float64
				for _, v := range q {
					sq += float64(v)
					if v < 30 {
						a++
					}
					if v < 20 {
						b++
					}
				}
				sq /= float64(n)
				s.sum[m] += sq
				s.sumSq[m] += sq * sq
				s.low30[m] = append(s.low30[m], a)
				s.low20[m] = append(s.low20[m], b)
			}
		}
		return s
	}
	c, ind := draw(true), draw(false)
	for _, x := range []struct {
		key  string
		c, i float64
	}{
		{"lowq30_count_correlation", correlation(c.low30[0], c.low30[1]), correlation(ind.low30[0], ind.low30[1])},
		{"lowq20_count_correlation", correlation(c.low20[0], c.low20[1]), correlation(ind.low20[0], ind.low20[1])},
	} {
		t.Logf("%s: coupled %.4f, independent %.4f, measured %.4f", x.key, x.c, x.i, want[x.key])
		if math.Abs(x.c-want[x.key]) > 0.03 {
			t.Errorf("%s: coupled mates %.4f, measured %.4f", x.key, x.c, want[x.key])
		}
		if math.Abs(x.i) > 0.02 {
			t.Errorf("%s: independent mates %.4f, want about 0", x.key, x.i)
		}
	}
	for m := range 2 {
		a, b := c.sum[m]/pairs, ind.sum[m]/pairs
		sd := math.Sqrt(ind.sumSq[m]/pairs - b*b) // of a read's mean Q
		tol := 4 * sd * math.Sqrt(2.0/pairs)
		t.Logf("read %d: mean Q coupled %.4f, independent %.4f (tolerance %.4f)", m+1, a, b, tol)
		if math.Abs(a-b) > tol {
			t.Errorf("read %d: mean Q %.4f with coupled classes, %.4f independent", m+1, a, b)
		}
	}
}

// correlation is the Pearson correlation of x and y.
func correlation(x, y []float64) float64 {
	var mx, my float64
	for i := range x {
		mx += x[i]
		my += y[i]
	}
	mx /= float64(len(x))
	my /= float64(len(y))
	var sxy, sxx, syy float64
	for i := range x {
		dx, dy := x[i]-mx, y[i]-my
		sxy += dx * dy
		sxx += dx * dx
		syy += dy * dy
	}
	return sxy / math.Sqrt(sxx*syy)
}
