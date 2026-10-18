//go:build !release && cgo

package genomics

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
)

// S3 output is spooled and finalized before upload: the uploaded objects
// carry the quality-score fingerprint (embedded, then gzip-compressed), and
// the manifest describes the uploaded bytes.
func TestS3_FingerprintVerifiesOnUploadedBytes(t *testing.T) {
	c := qualityFingerprintConfig(t.TempDir(), 60, 3, true)
	c.Seed = 17
	c.Workers = 3
	c.Output.Compression = "gzip"
	c.OutputDir = "s3://bk/fp/run"
	up, e := s3Run(t, c)

	m := e.Manifest()
	var want string
	fastqs := 0
	for _, o := range m.Outputs {
		data := up.objects["bk/"+o.Key]
		if sha256Hex(data) != o.SHA256 || int64(len(data)) != o.Bytes {
			t.Errorf("%s: manifest does not match the uploaded object", o.Key)
		}
		if o.Role != manifest.RoleFASTQ {
			continue
		}
		fastqs++
		if !strings.HasSuffix(o.Key, ".fastq.gz") {
			t.Errorf("object %s is not a .fastq.gz", o.Key)
		}
		gz := filepath.Join(t.TempDir(), "obj.fastq.gz")
		if err := os.WriteFile(gz, data, 0o600); err != nil {
			t.Fatal(err)
		}
		plain := filepath.Join(t.TempDir(), "obj.fastq")
		if err := os.WriteFile(plain, readMaybeGzip(t, gz), 0o600); err != nil {
			t.Fatal(err)
		}
		got, err := extractFileFingerprint(t, e, plain)
		if err != nil {
			t.Fatalf("%s: fingerprint does not verify on the uploaded bytes: %v", o.Key, err)
		}
		if want == "" {
			want = got
		} else if got != want {
			t.Errorf("objects carry different fingerprints: %q vs %q", got, want)
		}
	}
	if fastqs != 6 || len(up.objects) != len(m.Outputs)+1 {
		t.Errorf("%d FASTQ outputs, %d objects for %d outputs + manifest", fastqs, len(up.objects), len(m.Outputs))
	}
	if !m.Encrypted {
		t.Error("manifest does not record encrypted fingerprints")
	}
}
