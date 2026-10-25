//go:build !release

package illumina

import (
	"bytes"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// indelResults returns the directory of the indel measurements of p's run
// (scripts/acceptance/illumina_indels.sh).
func indelResults(p *Profile) string {
	sub := p.Run.Accession
	if len(p.Runs) > 0 {
		sub = "pooled"
	}
	return filepath.Join(resultsRoot, "indels", sub)
}

// binomialOK reports whether k successes in n trials are within z standard
// deviations of the probability p (with a floor of 3 events for tiny
// expectations).
func binomialOK(k, n int, p, z float64) bool {
	mean := float64(n) * p
	sd := math.Sqrt(mean * (1 - p))
	return math.Abs(float64(k)-mean) <= z*sd+3
}

// poissonCI returns the exact (Garwood) 95 % interval of a Poisson mean
// given k events.
func poissonCI(k int) (lo, hi float64) {
	cdf := func(k int, mu float64) float64 {
		if k < 0 {
			return 0
		}
		s := 0.0
		for i := 0; i <= k; i++ {
			lg, _ := math.Lgamma(float64(i + 1))
			s += math.Exp(float64(i)*math.Log(mu) - mu - lg)
		}
		return s
	}
	solve := func(f func(float64) bool, a, b float64) float64 {
		for range 200 {
			m := (a + b) / 2
			if f(m) {
				a = m
			} else {
				b = m
			}
		}
		return (a + b) / 2
	}
	top := float64(k) + 20*math.Sqrt(float64(k)+1) + 20
	if k > 0 {
		lo = solve(func(mu float64) bool { return 1-cdf(k-1, mu) < 0.025 }, 0, top)
	}
	hi = solve(func(mu float64) bool { return cdf(k, mu) > 0.025 }, 0, top)
	return lo, hi
}

// recoveryFactors returns the closed loop's errors generated per error
// found by kind (0 insertions, 1 deletions) and group (short, repeat, bg),
// from results/indels/indel_fit.tsv; 1 where it has none.
func recoveryFactors(t *testing.T) map[[2]string]float64 {
	t.Helper()
	defer func(d string) { resultsDir = d }(resultsDir)
	resultsDir = filepath.Join(resultsRoot, "indels")
	f := map[[2]string]float64{}
	for _, r := range readTSV(t, "indel_fit.tsv") {
		if r["section"] == "recovery" && r["parameter"] == "factor" {
			f[[2]string{r["key"], r["class"]}] = num(t, r["value"])
		}
	}
	return f
}

// The model's expected indel errors on the repeats and bases where a
// profile's run was measured (the overlaps of its mates,
// results/indels/<run>/indel_classes.tsv), less the share the measurement
// loses (the closed loop's recovery, by which fit_indels.py raised the
// rates), are inside the exact 95 % interval of the errors measured there,
// insertions and deletions apart. With TestIndelModel_ClassRates (the
// generator draws at the model's rates) and TestCalibration_IndelClosedLoop
// (generated reads measured the same way) this ties the generated rate to
// the measured one.
func TestCalibration_IndelRatesMatchNovaSeq(t *testing.T) {
	recovery := recoveryFactors(t)
	found := func(op, cls string, x int) float64 {
		g := "bg"
		switch {
		case cls == "hp" && x < LongRun:
			g = "short"
		case cls != "bg":
			g = "repeat"
		}
		if f, ok := recovery[[2]string{op, g}]; ok {
			return 1 / f
		}
		return 1
	}
	for _, p := range profiles {
		t.Run(p.Name, func(t *testing.T) {
			defer func(d string) { resultsDir = d }(resultsDir)
			resultsDir = indelResults(p)
			var obs, exp [2]float64
			var bases float64
			for _, r := range readTSV(t, "indel_classes.tsv") {
				op := 0
				if r["op"] == "D" {
					op = 1
				}
				n, k := num(t, r["opportunities"]), num(t, r["events"])
				var rate float64
				x, _ := strconv.Atoi(strings.TrimSuffix(r["bin"], "+"))
				switch r["class"] {
				case "hp":
					rate = rateAt([][]float64{p.Indels.HPIns, p.Indels.HPDel}[op], x-1)
				case "str":
					rate = rateAt([][]float64{p.Indels.STRIns, p.Indels.STRDel}[op], x-2)
				case "bg":
					rate = [2]float64{p.Indels.BgIns, p.Indels.BgDel}[op]
					bases = n
				}
				obs[op] += k
				exp[op] += n * rate * found(r["op"], r["class"], x)
			}
			for op, name := range []string{"insertions", "deletions"} {
				lo, hi := poissonCI(int(obs[op]))
				if exp[op] < lo || exp[op] > hi {
					t.Errorf("%s: model expects %.1f where %v were measured (95 %% interval %.1f-%.1f)", name, exp[op], obs[op], lo, hi)
				}
				t.Logf("%s: measured %v (%.3g per base, 95 %% %.3g-%.3g), model %.1f", name, obs[op], obs[op]/bases, lo/bases, hi/bases, exp[op])
			}
		})
	}
}

// The closed loop (scripts/acceptance/illumina_indels.sh with GENERATED):
// reads generated with a profile from GRCh38, aligned with minimap2 and
// measured like the real ones, have the profile's run's insertion and
// deletion rates: each measured generated rate's 95 % interval overlaps
// the real run's.
func TestCalibration_IndelClosedLoop(t *testing.T) {
	for _, p := range profiles {
		dir := filepath.Join(resultsRoot, "indels", "generated-"+p.Name)
		if _, err := os.Stat(filepath.Join(dir, "indel_summary.tsv")); err != nil {
			continue // not every profile was generated
		}
		t.Run(p.Name, func(t *testing.T) {
			defer func(d string) { resultsDir = d }(resultsDir)
			rates := func(d string) map[string][3]float64 {
				resultsDir = d
				out := map[string][3]float64{}
				for _, r := range readTSV(t, "indel_summary.tsv") {
					if r["ci_low"] != "" {
						out[r["key"]] = [3]float64{num(t, r["value"]), num(t, r["ci_low"]), num(t, r["ci_high"])}
					}
				}
				return out
			}
			gen, real := rates(dir), rates(indelResults(p))
			for _, k := range []string{"ins_rate_per_base", "del_rate_per_base"} {
				g, r := gen[k], real[k]
				if g[1] > r[2] || g[2] < r[1] {
					t.Errorf("%s: generated %.3g (%.3g-%.3g), real %.3g (%.3g-%.3g)", k, g[0], g[1], g[2], r[0], r[1], r[2])
				}
				t.Logf("%s: generated %.3g (%.3g-%.3g), real %.3g (%.3g-%.3g)", k, g[0], g[1], g[2], r[0], r[1], r[2])
			}
		})
	}
}

// testIndels is a model with rates high enough to test every tier.
var testIndels = IndelParams{
	HPIns:      []float64{0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1},
	HPDel:      []float64{0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.04, 0.05, 0.06, 0.07, 0.08, 0.12},
	STRIns:     []float64{0.004, 0.008, 0.012, 0.05, 0.1},
	STRDel:     []float64{0.003, 0.006, 0.009, 0.04, 0.08},
	BgIns:      0.001,
	BgDel:      0.0008,
	HPInsSize:  [2][]float64{{0.8, 0.2}, {0.5, 0.3, 0.2}},
	HPDelSize:  [2][]float64{{0.9, 0.1}, {0.6, 0.3, 0.1}},
	STRInsSize: []float64{0.7, 0.3},
	STRDelSize: []float64{0.75, 0.25},
	BgInsSize:  []float64{0.5, 0.3, 0.2},
	BgDelSize:  []float64{0, 0.6, 0.4},
	InsQual:    [4]float64{0, 0.3, 0.1, 0.6},
	DelQual:    [4]float64{0, 0.5, 0.2, 0.3},
}

// repeatTemplate returns a template of homopolymer runs of 1-20 bases and
// tandem repeats of 2-8 copies of 2-6 base units, separated by spacers
// that start no repeat, and the opportunities it holds by class.
func repeatTemplate(r *rand.Rand) []byte {
	var b bytes.Buffer
	spacer := func() { b.WriteString("ACGTCAGTGACT") }
	for h := 1; h <= 20; h++ {
		spacer()
		b.WriteString(strings.Repeat(string("ACGT"[r.IntN(4)]), h))
	}
	for c := 2; c <= 8; c++ {
		for _, u := range []string{"CA", "TTG", "AGAT", "AATGG", "AAGCTC"} {
			spacer()
			b.WriteString(strings.Repeat(u, c))
		}
	}
	spacer()
	return b.Bytes()
}

// apply reconstructs a read from its template and errors (with the
// inserted bases taken from the read), to check that the read is the
// template with exactly these errors.
func apply(t []byte, read []byte, edits []Indel, n int) []byte {
	var out []byte
	i, e := 0, 0
	for len(out) < n && (i < len(t) || e < len(edits)) {
		for e < len(edits) && edits[e].Pos == i && len(out) == edits[e].Cycle {
			ed := edits[e]
			e++
			if ed.Insertion {
				out = append(out, read[ed.Cycle:ed.Cycle+ed.Len]...)
				continue
			}
			i += ed.Len
		}
		if len(out) >= n || i >= len(t) {
			break
		}
		out = append(out, t[i])
		i++
	}
	return out
}

// The generator makes errors at each class's rate per opportunity: per
// homopolymer run by its length, per tandem repeat by its copies, per
// background base; thinned (short runs and tracts) and directly drawn
// alike; with the sizes of the model.
func TestIndelModel_ClassRates(t *testing.T) {
	m := NewIndelModel(testIndels)
	r := rand.New(rand.NewPCG(30, 1)) // #nosec G404 -- test stream
	tmpl := repeatTemplate(r)
	// opportunities per template, from the walk's own view
	type key struct {
		c   IndelClass
		rep int
		ins bool
	}
	events := map[key]int{}
	sizes := map[key][]int{}
	const reads = 40000
	for range reads {
		read, edits := m.Read(r, tmpl, -1, nil)
		_ = read
		for _, e := range edits {
			k := key{e.Class, e.Repeat, e.Insertion}
			events[k]++
			sz := e.Len
			if e.Class == IndelTandem {
				sz = e.Len / unitLen(tmpl, e.Pos)
			}
			s := sizes[k]
			for len(s) < sz {
				s = append(s, 0)
			}
			s[sz-1]++
			sizes[k] = s
		}
	}
	// opportunities: each run of length h and each tract of c copies once
	// per read (the template's runs and tracts)
	hpRuns, strTracts := map[int]int{}, map[int]int{}
	for i := 0; i < len(tmpl); i++ {
		if i == 0 || tmpl[i-1] != tmpl[i] {
			h := 1
			for i+h < len(tmpl) && tmpl[i+h] == tmpl[i] {
				h++
			}
			hpRuns[h]++
		}
		for k := strMinUnit; k <= strMaxUnit; k++ {
			if tractStart(tmpl, i, k) {
				strTracts[copiesAt(tmpl, i, k)]++
			}
		}
	}
	check := func(name string, k, n int, p float64) {
		t.Helper()
		if !binomialOK(k, n, p, 4) {
			t.Errorf("%s: %d errors in %d opportunities, want %.1f (p %.4g)", name, k, n, float64(n)*p, p)
		}
	}
	for h, runs := range hpRuns {
		n := runs * reads
		check("homopolymer "+strconv.Itoa(h)+" insertions", events[key{IndelHomopolymer, h, true}], n, rateAt(testIndels.HPIns, h-1))
		check("homopolymer "+strconv.Itoa(h)+" deletions", events[key{IndelHomopolymer, h, false}], n, rateAt(testIndels.HPDel, h-1))
	}
	for c, tr := range strTracts {
		n := tr * reads
		check("tandem "+strconv.Itoa(c)+" insertions", events[key{IndelTandem, c, true}], n, rateAt(testIndels.STRIns, c-2))
		check("tandem "+strconv.Itoa(c)+" deletions", events[key{IndelTandem, c, false}], n, rateAt(testIndels.STRDel, c-2))
	}
	// background: per template base (bases skipped by deletions are not
	// opportunities; their share is tiny here)
	n := len(tmpl) * reads
	check("background insertions", events[key{IndelBackground, 0, true}], n, testIndels.BgIns)
	check("background deletions", events[key{IndelBackground, 0, false}], n, testIndels.BgDel)
	// sizes of long-run slips and background errors
	dist := func(name string, got []int, want []float64) {
		t.Helper()
		tot := 0
		for _, x := range got {
			tot += x
		}
		for i, w := range want {
			k := 0
			if i < len(got) {
				k = got[i]
			}
			if !binomialOK(k, tot, w, 4) {
				t.Errorf("%s size %d: %d of %d, want %.3f", name, i+1, k, tot, w)
			}
		}
		if len(got) > len(want) {
			t.Errorf("%s: size %d drawn, the model has none above %d", name, len(got), len(want))
		}
	}
	var long []int
	for h := LongRun; h <= 20; h++ {
		for i, x := range sizes[key{IndelHomopolymer, h, false}] {
			for len(long) <= i {
				long = append(long, 0)
			}
			long[i] += x
		}
	}
	dist("long-run deletion", long, testIndels.HPDelSize[1])
	dist("background insertion", sizes[key{IndelBackground, 0, true}], testIndels.BgInsSize)
	dist("background deletion", sizes[key{IndelBackground, 0, false}], testIndels.BgDelSize)
}

// unitLen returns the unit length of the tandem repeat at t[i:].
func unitLen(t []byte, i int) int {
	for k := strMinUnit; k <= strMaxUnit; k++ {
		if tractStart(t, i, k) {
			return k
		}
	}
	return 1
}

// A read of n cycles is n bases, whatever its errors; it is its template
// with exactly the recorded errors (insertions of the run's base or unit in
// repeats); ReadOffset and TemplateUsed follow the errors; Qualities
// redraws the qualities at the errors from InsQual and DelQual.
func TestIndelModel_ReadsTemplateWithErrors(t *testing.T) {
	m := NewIndelModel(testIndels)
	r := rand.New(rand.NewPCG(30, 2)) // #nosec G404 -- test stream
	tmpl := repeatTemplate(r)
	var insQ, delQ [4]int
	withErrors := 0
	for range 20000 {
		const n = 151
		read, edits := m.Read(r, tmpl, n, nil)
		if len(read) != n {
			t.Fatalf("read of %d bases, want %d", len(read), n)
		}
		if got := apply(tmpl, read, edits, n); !bytes.Equal(got, read) {
			t.Fatalf("read is not its template with its errors %+v:\n%s\n%s", edits, read, got)
		}
		if used := TemplateUsed(n, edits); used > len(tmpl) {
			t.Fatalf("read used %d template bases of %d", used, len(tmpl))
		}
		for _, e := range edits {
			if e.Cycle < 0 || e.Cycle > n {
				t.Fatalf("error at cycle %d", e.Cycle)
			}
			if e.Insertion && e.Class == IndelHomopolymer {
				for _, b := range read[e.Cycle : e.Cycle+e.Len] {
					if b != tmpl[e.Pos] {
						t.Fatalf("homopolymer insertion %q into a run of %c", read[e.Cycle:e.Cycle+e.Len], tmpl[e.Pos])
					}
				}
			}
		}
		// The template base just past each error is where ReadOffset puts it.
		for _, e := range edits {
			x := e.Pos
			if !e.Insertion {
				x += e.Len
			}
			if deleted(edits, x) {
				continue
			}
			if off := ReadOffset(x, edits); off < n && x < len(tmpl) && read[off] != tmpl[x] {
				t.Fatalf("ReadOffset(%d) = %d: read has %c, template %c", x, off, read[off], tmpl[x])
			}
		}
		if len(edits) == 0 {
			continue
		}
		withErrors++
		qual := bytes.Repeat([]byte{QHigh}, n)
		m.Qualities(r, qual, edits)
		for _, e := range edits {
			if e.Insertion {
				for c := e.Cycle; c < e.Cycle+e.Len && c < n; c++ {
					insQ[binIndex(qual[c])]++
				}
			} else if e.Cycle < n && e.Cycle > 0 && !overlaps(edits, e) {
				delQ[binIndex(qual[e.Cycle])]++
			}
		}
	}
	if withErrors < 1000 {
		t.Fatalf("only %d reads with errors", withErrors)
	}
	for name, c := range map[string]struct {
		got  [4]int
		want [4]float64
	}{"inserted": {insQ, testIndels.InsQual}, "deletion flank": {delQ, testIndels.DelQual}} {
		tot := c.got[0] + c.got[1] + c.got[2] + c.got[3]
		for i := range c.want {
			if !binomialOK(c.got[i], tot, c.want[i], 4) {
				t.Errorf("%s bases: %d of %d at Q%d, want %.2f", name, c.got[i], tot, Bins[i], c.want[i])
			}
		}
	}
}

// deleted reports whether an error of edits deletes template base x.
func deleted(edits []Indel, x int) bool {
	for _, e := range edits {
		if !e.Insertion && e.Pos <= x && x < e.Pos+e.Len {
			return true
		}
	}
	return false
}

// overlaps reports whether another error of edits touches e's cycles.
func overlaps(edits []Indel, e Indel) bool {
	for _, o := range edits {
		if o != e && o.Cycle >= e.Cycle-1 && o.Cycle <= e.Cycle+e.Len {
			return true
		}
	}
	return false
}

// Deletions never leave a read short: with a template barely longer than
// the read the read is still n bases, made of template bases with errors;
// a read of the whole template (n < 0) changes length by its errors.
func TestIndelModel_ShortTemplates(t *testing.T) {
	m := NewIndelModel(testIndels)
	r := rand.New(rand.NewPCG(30, 3)) // #nosec G404 -- test stream
	tmpl := repeatTemplate(r)
	for i := range 5000 {
		n := 150 + i%3
		read, edits := m.Read(r, tmpl[:n+i%2], n, nil)
		if len(read) != n {
			t.Fatalf("read of %d bases from a template of %d, want %d", len(read), n+i%2, n)
		}
		if TemplateUsed(n, edits) > n+i%2 {
			t.Fatalf("read used %d bases of a %d-base template", TemplateUsed(n, edits), n+i%2)
		}
		whole, edits := m.Read(r, tmpl, -1, nil)
		if want := len(tmpl) - TemplateUsed(0, edits); len(whole) != want {
			t.Fatalf("whole-template read of %d bases, want %d", len(whole), want)
		}
	}
	// more extends the template when deletions need it
	asked := 0
	more := func(k int) []byte { asked += k; return bytes.Repeat([]byte("ACGTTGCA"), k/8+1)[:k] }
	for range 3000 {
		read, _ := m.Read(r, tmpl[:151], 151, more)
		if len(read) != 151 {
			t.Fatalf("read of %d bases with an extensible template", len(read))
		}
	}
	if asked == 0 {
		t.Error("the template was never extended")
	}
}

// Every profile has an indel model: rates in [0, 1), the background
// deletion sizes start at 2 (a one-base deletion is a homopolymer slip),
// distributions sum to 1, and the rates grow with the repeat (a long
// homopolymer slips far more often than a short one).
func TestIndelParams_Profiles(t *testing.T) {
	for _, p := range profiles {
		ip := p.Indels
		if !ip.Enabled() || NewIndelModel(ip) == nil {
			t.Fatalf("%s: no indel model", p.Name)
		}
		for name, tab := range map[string][]float64{"HPIns": ip.HPIns, "HPDel": ip.HPDel, "STRIns": ip.STRIns, "STRDel": ip.STRDel} {
			for i, x := range tab {
				if !(x >= 0 && x < 1) {
					t.Errorf("%s %s[%d] = %g", p.Name, name, i, x)
				}
			}
			if !(tab[len(tab)-1] > 100*tab[0]) {
				t.Errorf("%s %s: longest repeat %g, shortest %g; want a strong rise", p.Name, name, tab[len(tab)-1], tab[0])
			}
		}
		if ip.BgDelSize[0] != 0 {
			t.Errorf("%s: background one-base deletions %g", p.Name, ip.BgDelSize[0])
		}
		for name, d := range map[string][]float64{"HPInsSize0": ip.HPInsSize[0], "HPInsSize1": ip.HPInsSize[1], "HPDelSize0": ip.HPDelSize[0],
			"HPDelSize1": ip.HPDelSize[1], "STRInsSize": ip.STRInsSize, "STRDelSize": ip.STRDelSize, "BgInsSize": ip.BgInsSize, "BgDelSize": ip.BgDelSize,
			"InsQual": ip.InsQual[:], "DelQual": ip.DelQual[:]} {
			s := 0.0
			for _, x := range d {
				s += x
			}
			if math.Abs(s-1) > 1e-4 {
				t.Errorf("%s %s sums to %g", p.Name, name, s)
			}
		}
	}
}

// A model without errors makes none and draws nothing.
func TestIndelModel_Disabled(t *testing.T) {
	if m := NewIndelModel(IndelParams{}); m != nil {
		t.Fatal("empty parameters gave a model")
	}
	var m *IndelModel
	r := rand.New(rand.NewPCG(1, 1)) // #nosec G404 -- test stream
	read, edits := m.Read(r, []byte("ACGTACGT"), 5, nil)
	if string(read) != "ACGTA" || edits != nil {
		t.Fatalf("nil model read %q %v", read, edits)
	}
	if r.Uint64() != rand.New(rand.NewPCG(1, 1)).Uint64() { // #nosec G404 -- test stream
		t.Error("nil model drew from the stream")
	}
}

// Read is a pure function of the stream.
func TestIndelModel_Deterministic(t *testing.T) {
	m := NewIndelModel(profilePOOLED.Indels)
	tmpl := repeatTemplate(rand.New(rand.NewPCG(5, 5))) // #nosec G404 -- test data
	for s := range uint64(200) {
		a, ea := m.Read(rand.New(rand.NewPCG(s, 9)), tmpl, 151, nil) // #nosec G404 -- test stream
		b, eb := m.Read(rand.New(rand.NewPCG(s, 9)), tmpl, 151, nil) // #nosec G404 -- test stream
		if !bytes.Equal(a, b) || len(ea) != len(eb) {
			t.Fatalf("seed %d: reads differ", s)
		}
	}
}
