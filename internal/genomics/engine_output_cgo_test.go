//go:build !release && cgo

package genomics

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore"
)

// The quality-score fingerprint is embedded in the plain FASTQ and then the
// file is compressed: after decompressing the .fastq.gz the fingerprint must
// still verify, and the summary must still describe the reads.
func TestONT_GzipQualityFingerprintSurvivesCompression(t *testing.T) {
	dir := t.TempDir()
	c := qualityFingerprintConfig(dir, 200, 2, false)
	c.Platform = "ont"
	c.Seed = 99
	c.ReadLength = 400
	c.ONT.Channels = 32
	c.Output.Compression = "gzip"
	e := NewEngine(c)
	defer e.Close()
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	files := e.GeneratedFiles()
	if len(files) != 2 {
		t.Fatalf("GeneratedFiles = %v", files)
	}
	var want string
	for _, f := range files {
		if !strings.HasSuffix(f, ".fastq.gz") {
			t.Fatalf("%s is not a .fastq.gz", f)
		}
		if _, err := os.Stat(strings.TrimSuffix(f, ".gz")); err == nil {
			t.Errorf("plain %s left behind", filepath.Base(f))
		}
		plain := filepath.Join(t.TempDir(), "plain.fastq")
		if err := os.WriteFile(plain, readMaybeGzip(t, f), 0o600); err != nil {
			t.Fatal(err)
		}
		got, err := extractFileFingerprint(t, e, plain)
		if err != nil {
			t.Fatalf("%s: fingerprint does not verify after decompression: %v", filepath.Base(f), err)
		}
		if want == "" {
			want = got
		} else if got != want {
			t.Errorf("files carry different fingerprints: %q vs %q", got, want)
		}
	}
	if !strings.HasPrefix(want, "customer-7:dataset-42:") {
		t.Errorf("fingerprint %q", want)
	}

	// Headers, lengths, channels and times still agree with the summary,
	// which names the final .gz files. (Embedding perturbs some quality
	// scores after the summary is written, so mean qscores are compared
	// loosely here.)
	rows := readSummary(t, filepath.Join(dir, nanopore.SummaryFileName))
	if len(rows) != 2*c.ReadsPerBatch {
		t.Fatalf("%d summary rows", len(rows))
	}
	i := 0
	for _, f := range files {
		for _, rec := range parseFASTQ(t, string(readMaybeGzip(t, f))) {
			m := ontHeader.FindStringSubmatch(rec.header)
			if m == nil {
				t.Fatalf("header %q", rec.header)
			}
			r := rows[i]
			i++
			if r.id != m[1] || r.file != filepath.Base(f) || r.length != len(rec.seq) {
				t.Errorf("summary row %+v vs %s in %s", r, rec.header, filepath.Base(f))
			}
			if d := r.q - ontMeanQ(rec.qual); d > 3 || d < -3 {
				t.Errorf("read %s: mean_qscore %v far from recomputed %v", r.id, r.q, ontMeanQ(rec.qual))
			}
		}
	}
}
