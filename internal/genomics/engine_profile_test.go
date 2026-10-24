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

// Naming the default profile gives exactly the output of no profile; an
// unknown profile, or a profile with the linear model, is rejected before
// any output.
func TestIllumina_ProfileSelection(t *testing.T) {
	c := binnedConfig(t.TempDir())
	c.ReadsPerBatch = 300
	base := streamOutput(t, c, GenerateOpts{})
	c.Illumina.Profile = illumina.DefaultProfile
	if !bytes.Equal(streamOutput(t, c, GenerateOpts{}), base) {
		t.Error("the default profile named explicitly changed the output")
	}
	c.Illumina.Profile = strings.ToUpper(illumina.DefaultProfile)
	if !bytes.Equal(streamOutput(t, c, GenerateOpts{}), base) {
		t.Error("profile names are not case-insensitive")
	}
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
