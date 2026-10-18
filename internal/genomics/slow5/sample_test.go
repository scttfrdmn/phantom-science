//go:build !release

package slow5

import (
	"fmt"
	"os"
	"path/filepath"
	"testing"
)

// TestWriteSampleFiles writes a small synthetic flow cell for the acceptance
// script scripts/acceptance/slow5_validate.sh. It runs only when
// PHANTOM_SLOW5_SAMPLE_DIR is set, and writes into that directory:
//
//	chunks/<flowcell>_<runid>_<n>.blow5   zstd + svb-zd, 25 reads per file
//	codec_<rec>_<sig>.blow5               the same reads, one file per codec
//	reads.slow5                           the same reads as SLOW5 ASCII
//	read_ids.txt                          every read ID, one per line
func TestWriteSampleFiles(t *testing.T) {
	dir := os.Getenv("PHANTOM_SLOW5_SAMPLE_DIR")
	if dir == "" {
		t.Skip("set PHANTOM_SLOW5_SAMPLE_DIR to write sample files")
	}
	h := Header{
		ReadGroups: []ReadGroup{{
			"run_id":                  "8f1c2b3d4e5f60718293a4b5c6d7e8f901234567",
			"flow_cell_id":            "PAW12345",
			"flow_cell_product_code":  "FLO-PRO114M",
			"asic_id":                 "0004A30B00232BEC",
			"device_type":             "promethion",
			"exp_start_time":          "2024-03-01T10:00:00Z",
			"experiment_type":         "genomic_dna",
			"sample_frequency":        "5000",
			"sample_id":               "phantom_science",
			"sequencing_kit":          "sqk-lsk114",
			"protocol_group_id":       "phantom_science",
			"protocol_run_id":         "6a1f7c2e-4b9d-4e3a-8c5f-1d2e3f4a5b6c",
			"exp_script_purpose":      "sequencing_run",
			"experiment_duration_set": "4320",
			"hostname":                "phantom",
		}},
		Aux: ONTAuxFields,
	}
	recs := chunkRecords(100)
	if err := os.MkdirAll(dir, 0o750); err != nil {
		t.Fatal(err)
	}
	cw, err := NewChunkedWriter(filepath.Join(dir, "chunks"), h, ChunkOptions{ReadsPerFile: 25, BLOW5: BLOW5Options{RecordZstd, SignalSvbZd}, Overwrite: true})
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range recs {
		if err := cw.Write(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := cw.Close(); err != nil {
		t.Fatal(err)
	}
	for _, c := range allCodecs {
		name := "reads.slow5"
		if !c.ascii {
			name = fmt.Sprintf("codec_%s_%s.blow5", c.opts.RecordCompression, c.opts.SignalCompression)
		}
		if err := os.WriteFile(filepath.Join(dir, name), writeAll(t, c, h, recs), 0o600); err != nil {
			t.Fatal(err)
		}
	}
	var ids []byte
	for _, r := range recs {
		ids = append(ids, r.ReadID+"\n"...)
	}
	if err := os.WriteFile(filepath.Join(dir, "read_ids.txt"), ids, 0o600); err != nil {
		t.Fatal(err)
	}
	t.Logf("wrote %d chunk files and %d codec files to %s", len(cw.Files()), len(allCodecs), dir)
}
