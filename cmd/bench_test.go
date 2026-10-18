//go:build !release

package cmd

import (
	"bytes"
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/spf13/pflag"
)

// setBenchFlags sets bench flags for one test and restores them after.
func setBenchFlags(t *testing.T, kv map[string]string) {
	t.Helper()
	t.Cleanup(func() {
		benchCmd.Flags().VisitAll(func(f *pflag.Flag) {
			if f.Changed {
				_ = f.Value.Set(f.DefValue)
				f.Changed = false
			}
		})
		_ = benchCmd.Flags().Lookup("total-bases").Value.Set("")
	})
	for k, v := range kv {
		if err := benchCmd.Flags().Set(k, v); err != nil {
			t.Fatalf("set --%s: %v", k, err)
		}
	}
}

func TestBench_ReportsEveryWorkerCount(t *testing.T) {
	out := t.TempDir()
	setBenchFlags(t, map[string]string{
		"platform": "ont", "basecaller-tier": "fast", "total-bases": "150K",
		"workers": "1,3", "format": "json", "out": out,
	})
	var buf bytes.Buffer
	if err := runBench(benchCmd, &buf); err != nil {
		t.Fatal(err)
	}
	var r BenchReport
	if err := json.Unmarshal(buf.Bytes(), &r); err != nil {
		t.Fatalf("report is not JSON: %v\n%s", err, buf.String())
	}
	if len(r.Results) != 2 || r.Results[0].Workers != 1 || r.Results[1].Workers != 3 {
		t.Fatalf("results %+v", r.Results)
	}
	if !r.Identical || r.Scenario.Platform != "ont" || r.Scenario.Tier != "fast" || r.Scenario.TotalBases != 150000 {
		t.Errorf("report %+v", r)
	}
	for _, x := range r.Results {
		if x.Bases < 150000 || x.Reads <= 0 || x.Bytes <= x.Bases || x.Seconds <= 0 ||
			x.ReadsPerSecond <= 0 || x.BasesPerSecond <= 0 || x.MBPerSecond <= 0 {
			t.Errorf("implausible result %+v", x)
		}
		// The kept run is a complete, verifiable output.
		m, err := manifest.Read(filepath.Join(x.Dir, manifest.FileName))
		if err != nil {
			t.Fatal(err)
		}
		if err := m.VerifyDir(x.Dir); err != nil {
			t.Error(err)
		}
		if m.Totals.Bases != x.Bases {
			t.Errorf("manifest bases %d, report %d", m.Totals.Bases, x.Bases)
		}
	}
	saved, err := os.ReadFile(filepath.Join(out, "bench.json"))
	if err != nil || !bytes.Equal(saved, buf.Bytes()) {
		t.Errorf("bench.json not written with the report: %v", err)
	}
}

func TestBench_TableAndErrors(t *testing.T) {
	setBenchFlags(t, map[string]string{"total-bases": "30K", "workers": "2"})
	var buf bytes.Buffer
	if err := runBench(benchCmd, &buf); err != nil {
		t.Fatal(err)
	}
	s := buf.String()
	for _, want := range []string{"platform=illumina", "reads/s", "bases/s", "MB/s", "identical across worker counts: yes"} {
		if !strings.Contains(s, want) {
			t.Errorf("table lacks %q:\n%s", want, s)
		}
	}
	for flag, val := range map[string]string{"workers": "0,2", "format": "xml"} {
		setBenchFlags(t, map[string]string{flag: val})
		if err := runBench(benchCmd, &bytes.Buffer{}); err == nil {
			t.Errorf("--%s %s accepted", flag, val)
		}
	}
}
