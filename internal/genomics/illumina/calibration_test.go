//go:build !release

package illumina

import (
	"bufio"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// These tests compare the model with the measurements of real NovaSeq 6000
// reads (ENA DRR276626, 200,000 pairs of 2×151) committed by
// scripts/acceptance/illumina_calibrate.sh. The raw reads are not needed.

var resultsDir = filepath.Join("..", "..", "..", "scripts", "acceptance", "calibrate", "illumina", "results")

// readTSV returns the rows of a results table as maps from column name.
func readTSV(t *testing.T, name string) []map[string]string {
	t.Helper()
	f, err := os.Open(filepath.Join(resultsDir, name))
	if err != nil {
		t.Fatalf("measured table: %v", err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Scan()
	head := strings.Split(sc.Text(), "\t")
	var rows []map[string]string
	for sc.Scan() {
		row := map[string]string{}
		for i, v := range strings.Split(sc.Text(), "\t") {
			row[head[i]] = v
		}
		rows = append(rows, row)
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return rows
}

func num(t *testing.T, s string) float64 {
	t.Helper()
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		t.Fatalf("%q: %v", s, err)
	}
	return v
}

// realMate holds the measured statistics of one mate.
type realMate struct {
	perCycle   [][4]float64 // bin fractions per cycle
	lowCDF     []float64    // P(reads with <= k bases below Q37)
	same       [4]float64   // P(bin at c+1 = bin at c | bin)
	meanQ, q30 float64
}

func loadReal(t *testing.T) [2]realMate {
	var m [2]realMate
	for _, r := range readTSV(t, "cycle_quals.tsv") {
		i := int(num(t, r["mate"])) - 1
		n := num(t, r["bases"])
		var f [4]float64
		for b, q := range Bins {
			f[b] = num(t, r["Q"+strconv.Itoa(int(q))]) / n
		}
		m[i].perCycle = append(m[i].perCycle, f)
	}
	counts := [2][]float64{make([]float64, 152), make([]float64, 152)}
	for _, r := range readTSV(t, "read_composition.tsv") {
		i := int(num(t, r["mate"])) - 1
		k := int(num(t, r["n_Q2"]) + num(t, r["n_Q11"]) + num(t, r["n_Q25"]))
		counts[i][k] += num(t, r["reads"])
	}
	for i := range m {
		total, acc := 0.0, 0.0
		for _, v := range counts[i] {
			total += v
		}
		for _, v := range counts[i] {
			acc += v
			m[i].lowCDF = append(m[i].lowCDF, acc/total)
		}
	}
	var from, same [2][4]float64
	for _, r := range readTSV(t, "transitions.tsv") {
		i := int(num(t, r["mate"])) - 1
		a, b := binIndex(byte(num(t, r["from"]))), binIndex(byte(num(t, r["to"])))
		from[i][a] += num(t, r["count"])
		if a == b {
			same[i][a] += num(t, r["count"])
		}
	}
	for i := range m {
		for b := range 4 {
			if from[i][b] > 0 {
				m[i].same[b] = same[i][b] / from[i][b]
			}
		}
	}
	for _, r := range readTSV(t, "fastq_summary.tsv") {
		i := int(num(t, r["mate"])) - 1
		switch r["key"] {
		case "mean_q":
			m[i].meanQ = num(t, r["value"])
		case "frac_ge_q30":
			m[i].q30 = num(t, r["value"])
		case "alphabet_phred":
			for _, s := range strings.Split(r["value"], ",") {
				if binIndex(byte(num(t, s))) < 0 {
					t.Errorf("read %d: measured quality %s is not in Bins %v", i+1, s, Bins)
				}
			}
		}
	}
	return m
}

// Generated reads have the measured per-cycle bin frequencies (mean |Δ|
// over cycles < 0.0025 per bin, every cycle within 0.015: the smooth curve
// does not follow single-cycle features such as read 2 cycles 148-149,
// about 0.01 worse than their neighbours), the measured
// distribution of low-quality bases per read (CDF within 0.02), the
// measured adjacent-cycle repeat rates of Q11 and Q25 (within 0.01), mean
// quality (within 0.05) and ≥ Q30 fraction (within 0.003); read 2 is worse
// than read 1, as measured.
func TestCalibration_QualitiesMatchNovaSeq(t *testing.T) {
	real := loadReal(t)
	loadReadCells(t)
	const reads = 100000
	var genMean [2]float64
	for i, p := range []QualityParams{Read1, Read2} {
		rm := real[i]
		n := len(rm.perCycle)
		if n != 151 {
			t.Fatalf("read %d: %d measured cycles, want 151", i+1, n)
		}
		s := generate(p, uint64(10+i), reads, n, func(r *rand.Rand, n int) []byte { return humanBasesMate(r, n, i) })
		for b := 1; b < 4; b++ {
			sum, worst, at := 0.0, 0.0, 0
			for c := range n {
				d := math.Abs(float64(s.perCycle[c][b])/reads - rm.perCycle[c][b])
				sum += d
				if d > worst {
					worst, at = d, c+1
				}
			}
			mean := sum / float64(n)
			t.Logf("read %d Q%d: mean |Δ| %.5f, max %.5f at cycle %d", i+1, Bins[b], mean, worst, at)
			if mean > 0.0025 || worst > 0.015 {
				t.Errorf("read %d Q%d per-cycle fractions: mean |Δ| %.5f, max %.5f at cycle %d", i+1, Bins[b], mean, worst, at)
			}
		}
		acc, worst, at := 0.0, 0.0, 0
		for k := 0; k <= n; k++ {
			acc += float64(s.lowCounts[k]) / reads
			if d := math.Abs(acc - rm.lowCDF[k]); d > worst {
				worst, at = d, k
			}
		}
		t.Logf("read %d: low-quality bases per read, max |ΔCDF| %.4f at %d", i+1, worst, at)
		if worst > 0.02 {
			t.Errorf("read %d: per-read low-quality CDF differs by %.4f at %d bases", i+1, worst, at)
		}
		for b := 1; b <= 2; b++ {
			got := float64(s.same[b]) / float64(s.from[b])
			t.Logf("read %d: P(Q%d|Q%d) %.4f, measured %.4f", i+1, Bins[b], Bins[b], got, rm.same[b])
			if math.Abs(got-rm.same[b]) > 0.01 {
				t.Errorf("read %d: P(Q%d|Q%d) = %.4f, measured %.4f", i+1, Bins[b], Bins[b], got, rm.same[b])
			}
		}
		bases := float64(reads * n)
		genMean[i] = s.sum / bases
		q30 := 0.0
		for c := range n {
			q30 += float64(s.perCycle[c][3])
		}
		q30 /= bases
		t.Logf("read %d: mean Q %.3f (measured %.3f), >=Q30 %.4f (measured %.4f)", i+1, genMean[i], rm.meanQ, q30, rm.q30)
		if math.Abs(genMean[i]-rm.meanQ) > 0.05 || math.Abs(q30-rm.q30) > 0.003 {
			t.Errorf("read %d: mean Q %.3f, >=Q30 %.4f; measured %.3f, %.4f", i+1, genMean[i], q30, rm.meanQ, rm.q30)
		}
	}
	if !(real[0].meanQ > real[1].meanQ) || !(genMean[0] > genMean[1]) {
		t.Errorf("mean Q: measured R1 %.3f R2 %.3f, model R1 %.3f R2 %.3f; want R1 > R2 in both",
			real[0].meanQ, real[1].meanQ, genMean[0], genMean[1])
	}
}

// The substitution model reproduces the measured rates: per bin and block
// of 30 cycles, the mismatch rate against GRCh38 less the sample's variant
// rate is within 15 % of Rate (and within 5 % over all cycles); the Q37 rate
// is the mate-overlap rate. The final cycle (151) is logged, not checked:
// its Q25 rate is about 1.3 times the curve's, which the model does not
// reproduce.
func TestCalibration_ErrorRatesMatchNovaSeq(t *testing.T) {
	variant := math.NaN()
	for _, r := range readTSV(t, "fit.tsv") {
		if r["param"] == "variant_rate" {
			variant = num(t, r["value"])
		}
	}
	if math.IsNaN(variant) {
		t.Fatal("fit.tsv has no variant_rate")
	}
	type acc struct{ aligned, mism, model float64 }
	blocks := map[[2]int]*acc{}
	all := map[int]*acc{}
	for _, r := range readTSV(t, "errors_cycle.tsv") {
		q, c := int(num(t, r["phred"])), int(num(t, r["cycle"]))
		if q != QLow && q != QMid {
			continue
		}
		a, k := num(t, r["aligned"]), num(t, r["mismatches"])
		m := NovaSeqErrors.Rate(byte(q), c-1) * a
		for _, x := range []*acc{getAcc(blocks, [2]int{q, min(c-1, 150) / 30}), getAcc(all, q)} {
			x.aligned += a
			x.mism += k
			x.model += m
		}
	}
	check := func(label string, x *acc, tol float64) {
		got := x.mism/x.aligned - variant
		want := x.model / x.aligned
		t.Logf("%s: measured %.5f (n=%.0f), model %.5f", label, got, x.aligned, want)
		if math.Abs(got-want) > tol*want {
			t.Errorf("%s: measured substitution rate %.5f, model %.5f", label, got, want)
		}
	}
	for q := range all {
		check("Q"+strconv.Itoa(q)+" all cycles", all[q], 0.05)
		for b := range 5 {
			check("Q"+strconv.Itoa(q)+" cycles from "+strconv.Itoa(30*b+1), getAcc(blocks, [2]int{q, b}), 0.15)
		}
		x := getAcc(blocks, [2]int{q, 5})
		t.Logf("Q%d final cycle: measured %.5f, model %.5f (not checked)", q, x.mism/x.aligned-variant, x.model/x.aligned)
	}
	for _, r := range readTSV(t, "errors_quality.tsv") {
		if r["mate"] == "both" && r["phred"] == "37" {
			if got := num(t, r["overlap_rate"]); math.Abs(got-NovaSeqErrors.High) > 1e-6 {
				t.Errorf("Q37 rate %v, measured overlap rate %v", NovaSeqErrors.High, got)
			}
		}
	}
}

func getAcc[K comparable, V any](m map[K]*V, k K) *V {
	if m[k] == nil {
		m[k] = new(V)
	}
	return m[k]
}
