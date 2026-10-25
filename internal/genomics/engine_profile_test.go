//go:build !release

package genomics

// Tests of the Illumina run profiles (illumina.profile, #35).

import (
	"bytes"
	"math"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/illumina"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// The profile and the insert size are resolved by the rule documented in
// package illumina (profiles.go): no profile gives the default profile
// (novaseq6000-pooled) with lognormal 455 ± 100 inserts; a profile set
// explicitly, the default included, brings its run's inserts; a set
// insert_size always wins (its empty fields defaulted).
func TestIllumina_ProfileResolution(t *testing.T) {
	ln := func(m, sd float64) types.InsertSizeSpec {
		return types.InsertSizeSpec{Dist: types.InsertSizeLognormal, Mean: m, SD: sd}
	}
	normal300 := types.InsertSizeSpec{Dist: types.InsertSizeNormal, Mean: 300, SD: 30}
	for _, tc := range []struct {
		name    string
		spec    types.IlluminaSpec
		profile string
		insert  types.InsertSizeSpec
	}{
		{"unset", types.IlluminaSpec{}, "novaseq6000-pooled", ln(455, 100)},
		{"explicit pooled", types.IlluminaSpec{Profile: "novaseq6000-pooled"}, "novaseq6000-pooled", ln(465, 140)},
		{"explicit drr394916", types.IlluminaSpec{Profile: "novaseq6000-drr394916"}, "novaseq6000-drr394916", ln(560, 115)},
		{"explicit drr394916, other case", types.IlluminaSpec{Profile: "NovaSeq6000-DRR394916"}, "novaseq6000-drr394916", ln(560, 115)},
		{"explicit drr276626", types.IlluminaSpec{Profile: "novaseq6000-drr276626"}, "novaseq6000-drr276626", ln(455, 100)},
		{"explicit err15600747", types.IlluminaSpec{Profile: "novaseq6000-err15600747"}, "novaseq6000-err15600747", ln(390, 150)},
		{"insert_size, no profile", types.IlluminaSpec{InsertSize: types.InsertSizeSpec{Mean: 300, SD: 30}},
			"novaseq6000-pooled", normal300},
		{"insert_size wins over drr394916", types.IlluminaSpec{Profile: "novaseq6000-drr394916", InsertSize: types.InsertSizeSpec{Mean: 300, SD: 30}},
			"novaseq6000-drr394916", normal300},
		{"insert_size wins over err15600747", types.IlluminaSpec{Profile: "novaseq6000-err15600747", InsertSize: ln(250, 40)},
			"novaseq6000-err15600747", ln(250, 40)},
		{"partial insert_size: defaults, not the profile's", types.IlluminaSpec{Profile: "novaseq6000-drr394916", InsertSize: types.InsertSizeSpec{Max: 900}},
			"novaseq6000-drr394916", types.InsertSizeSpec{Dist: types.InsertSizeLognormal, Mean: 455, SD: 100, Max: 900}},
		{"linear model", types.IlluminaSpec{QualityModel: types.IlluminaQualityLinear}, "novaseq6000-pooled", ln(455, 100)},
	} {
		t.Run(tc.name, func(t *testing.T) {
			c := binnedConfig(t.TempDir())
			c.Illumina = tc.spec
			e := NewEngine(c)
			defer e.Close()
			if err := e.checkIllumina(); err != nil {
				t.Fatal(err)
			}
			p, err := e.illuminaProfile()
			if err != nil || p.Name != tc.profile {
				t.Errorf("profile %v (err %v), want %s", p, err, tc.profile)
			}
			got, err := e.insertSizes()
			if err != nil || got != tc.insert {
				t.Errorf("insert size %+v (err %v), want %+v", got, err, tc.insert)
			}
		})
	}
}

// Naming the default profile gives the qualities, errors and mate coupling
// of no profile and differs only in the inserts: single-end output is
// identical, and paired output is identical once the inserts are set to the
// unset default (lognormal 455 ± 100) but not without. An unknown profile,
// or a profile with the linear model, is rejected before any output.
func TestIllumina_ProfileSelection(t *testing.T) {
	for _, paired := range []bool{false, true} {
		c := binnedConfig(t.TempDir())
		c.ReadsPerBatch, c.PairedEnd = 300, paired
		base := streamOutput(t, c, GenerateOpts{})
		for _, name := range []string{illumina.DefaultProfile, strings.ToUpper(illumina.DefaultProfile)} {
			c.Illumina = types.IlluminaSpec{Profile: name}
			if same := bytes.Equal(streamOutput(t, c, GenerateOpts{}), base); same == paired {
				t.Errorf("paired %v, profile %s named: output equal to no profile's %v, want %v", paired, name, same, !paired)
			}
			c.Illumina.InsertSize = types.InsertSizeSpec{Dist: types.InsertSizeLognormal, Mean: types.DefaultInsertSizeMean, SD: types.DefaultInsertSizeSD}
			if !bytes.Equal(streamOutput(t, c, GenerateOpts{}), base) {
				t.Errorf("paired %v, profile %s with the default inserts: output differs from no profile's", paired, name)
			}
		}
	}
	c := binnedConfig(t.TempDir())
	c.ReadsPerBatch = 300
	for _, tc := range []struct {
		spec types.IlluminaSpec
		want string
	}{
		{types.IlluminaSpec{Profile: "novaseq6000-nosuchrun"}, "unknown illumina profile"},
		{types.IlluminaSpec{Profile: illumina.DefaultProfile, QualityModel: types.IlluminaQualityLinear}, "needs quality_model binned"},
	} {
		c.Illumina = tc.spec
		e := NewEngine(c)
		err := e.GenerateToWriter(&bytes.Buffer{}, GenerateOpts{})
		e.Close()
		if err == nil || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("%+v: err = %v, want %q", tc.spec, err, tc.want)
		}
	}
}

// Each profile gives its run's quality level (mean Q of each mate within
// 0.1 of the run's, measured on 4,000 pairs of uniform bases) and, unless
// illumina.insert_size is set, its run's insert sizes.
func TestIllumina_ProfilesFollowTheirRuns(t *testing.T) {
	for _, p := range illumina.Profiles() {
		t.Run(p.Name, func(t *testing.T) {
			c := binnedConfig(t.TempDir())
			c.ReadLength, c.ReadsPerBatch = p.Run.ReadLength, 4000
			c.Illumina = types.IlluminaSpec{Profile: p.Name, Errors: types.IlluminaErrorsNone}
			e := NewEngine(c)
			defer e.Close()
			s := &recordingSampler{}
			e.fragments = s.sample
			var buf bytes.Buffer
			if err := e.GenerateToWriter(&buf, GenerateOpts{}); err != nil {
				t.Fatal(err)
			}
			recs := parseFASTQ(t, buf.String())
			var q [2]float64
			for i, r := range recs {
				for j := range len(r.qual) {
					q[i%2] += float64(r.qual[j] - '!')
				}
			}
			for m := range 2 {
				got := q[m] / float64(len(recs)/2*c.ReadLength)
				t.Logf("read %d: mean Q %.3f, run %.3f", m+1, got, p.Stats.MeanQ[m])
				if math.Abs(got-p.Stats.MeanQ[m]) > 0.1 {
					t.Errorf("read %d: mean Q %.3f, the run's %.3f", m+1, got, p.Stats.MeanQ[m])
				}
			}
			var sum float64
			for _, f := range s.frags {
				sum += float64(len(f))
			}
			mean := sum / float64(len(s.frags))
			t.Logf("insert mean %.1f, profile %v (run %.1f)", mean, p.InsertSize, p.Stats.InsertMean)
			if math.Abs(mean-p.InsertSize.Mean) > 4*p.InsertSize.SD/math.Sqrt(float64(len(s.frags)))+1 {
				t.Errorf("insert mean %.1f, profile's %v", mean, p.InsertSize.Mean)
			}

			c.Illumina.InsertSize = types.InsertSizeSpec{Mean: 300, SD: 30}
			e2 := NewEngine(c)
			defer e2.Close()
			s2 := &recordingSampler{}
			e2.fragments = s2.sample
			if err := e2.GenerateToWriter(&bytes.Buffer{}, GenerateOpts{}); err != nil {
				t.Fatal(err)
			}
			sum = 0
			for _, f := range s2.frags {
				sum += float64(len(f))
			}
			if m := sum / float64(len(s2.frags)); math.Abs(m-300) > 3 {
				t.Errorf("insert_size set: mean %.1f, want 300", m)
			}
		})
	}
}
