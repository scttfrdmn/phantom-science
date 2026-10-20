//go:build !release

package genomics

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Regenerate reproduces every read of every file, in order, from the
// manifest's file list or from the configured plan, and stops at limits.
func TestRegenerate_ReproducesRun(t *testing.T) {
	for name, mk := range map[string]func(string) types.SequencerConfig{
		"illumina": func(dir string) types.SequencerConfig {
			c := binnedConfig(dir)
			c.ReadsPerBatch, c.TotalReads = 70, 150
			return c
		},
		"ont_bases": func(dir string) types.SequencerConfig {
			c := scaleONTConfig(dir)
			c.TotalBatches, c.TotalBases, c.FlowCells = 0, 30000, 2
			return c
		},
	} {
		t.Run(name, func(t *testing.T) {
			dir := t.TempDir()
			c := mk(dir)
			tree, e := runTree(t, c)
			var want []string
			m := e.Manifest()
			var files []RunFile
			unit := map[int]int{}
			mates := len(e.mates())
			for _, o := range m.Outputs {
				if o.Role != manifest.RoleFASTQ {
					continue
				}
				for _, r := range parseFASTQ(t, string(readTreeFile(t, tree, o.Path))) {
					want = append(want, r.seq)
				}
				n := unit[o.FlowCell]
				unit[o.FlowCell]++
				files = append(files, RunFile{FlowCell: o.FlowCell, Batch: n/mates + 1, Mate: n%mates + 1, Reads: int(o.Reads)})
			}
			regen := func(opts RegenOptions) ([]string, int) {
				var got []string
				r := NewEngineLog(mk(""), nil)
				n, err := r.Regenerate(opts, func(_ RunFile, seq []byte) error {
					got = append(got, string(seq))
					return nil
				})
				if err != nil {
					t.Fatal(err)
				}
				return got, n
			}
			got, n := regen(RegenOptions{Files: files})
			if n != len(files) || strings.Join(got, "\n") != strings.Join(want, "\n") {
				t.Errorf("manifest files: %d files, %d reads; want %d files, %d reads, identical", n, len(got), len(files), len(want))
			}
			// The plan regenerates whole units, a superset of a run cut by
			// its bases target.
			planned, _ := regen(RegenOptions{})
			have := map[string]bool{}
			for _, seq := range planned {
				have[seq] = true
			}
			for i, seq := range want {
				if !have[seq] {
					t.Fatalf("planned regeneration (%d reads) lacks read %d of the run", len(planned), i)
				}
			}
			if got, n := regen(RegenOptions{Files: files, MaxReads: 5}); len(got) != 5 || n != 1 {
				t.Errorf("MaxReads 5: %d reads from %d files", len(got), n)
			}
			if _, n := regen(RegenOptions{Files: files, MaxFiles: 1}); n != 1 {
				t.Errorf("MaxFiles 1: %d files", n)
			}
		})
	}
}

func readTreeFile(t *testing.T, tree map[string][]byte, rel string) []byte {
	t.Helper()
	data := tree[rel]
	if strings.HasSuffix(rel, ".gz") {
		p := filepath.Join(t.TempDir(), "f.gz")
		if err := os.WriteFile(p, data, 0o600); err != nil {
			t.Fatal(err)
		}
		return readMaybeGzip(t, p)
	}
	return data
}
