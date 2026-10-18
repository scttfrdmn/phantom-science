//go:build !release

package slow5

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// TestSlow5toolsInterop checks our files against the reference tools. It
// runs only when SLOW5TOOLS names a slow5tools binary built with zstd
// support (make zstd=1), e.g.
//
//	SLOW5TOOLS=$(command -v slow5tools) go test -tags '!release' -run Slow5tools ./internal/genomics/slow5/
//
// For every codec it checks that slow5tools accepts the file (quickcheck,
// stats, index, get), that slow5tools' ASCII view of it is byte-identical to
// our SLOW5 writer's output, and that files slow5tools re-encodes from ours
// decode back to the same records here.
func TestSlow5toolsInterop(t *testing.T) {
	tool := os.Getenv("SLOW5TOOLS")
	if tool == "" {
		t.Skip("set SLOW5TOOLS to a slow5tools binary to run the interoperability test")
	}
	run := func(args ...string) []byte {
		t.Helper()
		var stdout, stderr bytes.Buffer
		cmd := exec.Command(tool, args...) // #nosec G204 -- test-only, path from the developer's environment
		cmd.Stdout, cmd.Stderr = &stdout, &stderr
		if err := cmd.Run(); err != nil {
			t.Fatalf("slow5tools %s: %v\n%s", strings.Join(args, " "), err, stderr.String())
		}
		return stdout.Bytes()
	}

	h := multiGroupHeader()
	r := testRand(11)
	var recs []*Record
	for i, n := range []int{0, 1, 2, 5, 100, 1000, 20000, 3} {
		rec := fullRecord(r, i, uint32(i%3), n)
		// Values whose %f form is exact, so slow5tools' ASCII rendering
		// equals ours byte for byte.
		rec.Offset, rec.Range = -243, 748.5
		rec.Aux["median_before"] = 210.125
		rec.Aux["a_f64"] = []float64{1.5, -2.25}
		recs = append(recs, rec)
	}
	recs[7].Aux = map[string]any{"channel_number": "12"} // mostly missing
	dir := t.TempDir()
	ascii := writeAll(t, allCodecs[0], h, recs)

	for _, c := range allCodecs {
		ext := ".blow5"
		if c.ascii {
			ext = ".slow5"
		}
		path := filepath.Join(dir, c.name+ext)
		if err := os.WriteFile(path, writeAll(t, c, h, recs), 0o600); err != nil {
			t.Fatal(err)
		}
		run("quickcheck", path)
		stats := string(run("stats", path))
		if !strings.Contains(stats, fmt.Sprintf("number of records\t%d", len(recs))) {
			t.Errorf("%s: stats does not report %d records:\n%s", c.name, len(recs), stats)
		}
		run("index", path)
		got := run("get", path, recs[4].ReadID, recs[0].ReadID, "--to", "slow5")
		if lines := strings.Split(strings.TrimSpace(string(got)), "\n"); !strings.HasPrefix(lines[len(lines)-2], recs[4].ReadID+"\t") ||
			!strings.HasPrefix(lines[len(lines)-1], recs[0].ReadID+"\t") {
			t.Errorf("%s: get returned the wrong reads", c.name)
		}
		view := run("view", path, "--to", "slow5")
		if !bytes.Equal(view, ascii) {
			t.Errorf("%s: slow5tools view differs from our SLOW5 at byte %d", c.name, firstDiff(view, ascii))
		}
		// slow5tools re-encodes; we must read its output. slow5tools 1.5.0
		// cannot svb-zd-compress a zero-length signal parsed from SLOW5
		// ASCII ("slow5_ptr_compress_solo: Argument 'ptr' cannot be NULL",
		// its parser leaves raw_signal NULL), so ASCII input is re-encoded
		// without signal compression.
		reencode := [][2]string{{"zstd", "svb-zd"}, {"zlib", "svb-zd"}, {"none", "none"}}
		if c.ascii {
			reencode = [][2]string{{"zstd", "none"}, {"none", "none"}}
		}
		for _, re := range reencode {
			out := filepath.Join(dir, c.name+"_"+re[0]+"_"+re[1]+".blow5")
			run("view", path, "-c", re[0], "-s", re[1], "-o", out)
			b, err := os.ReadFile(out)
			if err != nil {
				t.Fatal(err)
			}
			gotH, back := readAllBytes(t, b)
			if d := headerDiff(h, gotH); d != "" {
				t.Errorf("%s -> %v: header %s", c.name, re, d)
			}
			for i := range recs {
				if d := recordDiff(recs[i], back[i]); d != "" {
					t.Errorf("%s -> %v record %d: %s", c.name, re, i, d)
				}
			}
		}
	}

	// A chunked run: every file passes quickcheck and the merged view has
	// every read.
	chunkDir := filepath.Join(dir, "chunks")
	cw := writeChunks(t, chunkDir, ChunkOptions{ReadsPerFile: 3, BLOW5: BLOW5Options{RecordZstd, SignalSvbZd}}, chunkRecords(10))
	for _, f := range cw.Files() {
		run("quickcheck", f)
	}
	merged := filepath.Join(dir, "merged.blow5")
	run("merge", chunkDir, "-o", merged)
	if _, back := readAllBytes(t, mustRead(t, merged)); len(back) != 10 {
		t.Errorf("merged chunks have %d reads, want 10", len(back))
	}
}

func mustRead(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return b
}
