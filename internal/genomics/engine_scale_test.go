//go:build !release

package genomics

// Output-size targets, flow cells, parallel generation, the run manifest
// and S3 output (#16).

import (
	"bytes"
	"context"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strings"
	"sync"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/nanopore"
	"github.com/scttfrdmn/phantom-science/internal/storage/s3sink"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// runTree runs Generate and returns every file under the output directory
// (relative path → content), failing if a spool directory is left behind.
func runTree(t *testing.T, c types.SequencerConfig) (map[string][]byte, *Engine) {
	t.Helper()
	e := NewEngine(c)
	t.Cleanup(e.Close)
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	tree := map[string][]byte{}
	err := filepath.WalkDir(c.OutputDir, func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if strings.HasPrefix(d.Name(), ".phantom-spool-") || strings.HasSuffix(d.Name(), ".part") {
			t.Errorf("spool left behind: %s", path)
		}
		if d.IsDir() {
			return nil
		}
		data, err := os.ReadFile(path)
		if err != nil {
			return err
		}
		rel, _ := filepath.Rel(c.OutputDir, path)
		tree[filepath.ToSlash(rel)] = data
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	return tree, e
}

// trainedModel saves a small trained adaptive HMM model and returns its path.
func trainedModel(t *testing.T) string {
	t.Helper()
	model := authenticity.NewAdaptiveHMMModel(3)
	if err := model.TrainFromSequence(makeEngineTestSequence(5000)); err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(t.TempDir(), "m.hmm.json")
	if err := model.Save(path); err != nil {
		t.Fatal(err)
	}
	return path
}

func scaleONTConfig(dir string) types.SequencerConfig {
	c := ontConfig(dir)
	c.ReadsPerBatch = 12
	c.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthLognormal, N50: 900, Median: 600}
	c.ErrorModel.Tier = "hac"
	c.Output.Compression = "gzip"
	return c
}

// Output does not depend on the number of workers: files, the sequencing
// summary and the manifest are byte-identical for 1 and 8 workers, also
// when the bases target cuts a batch and with an HMM model loaded.
func TestScale_WorkersByteIdentical(t *testing.T) {
	model := trainedModel(t)
	cases := map[string]func(*types.SequencerConfig){
		"illumina fixed paired": func(c *types.SequencerConfig) { c.PairedEnd = true; c.TotalBatches = 9 },
		// The #27 default: binned qualities, quality-matched errors and
		// paired fragments regenerated in each mate file.
		"illumina binned paired": func(c *types.SequencerConfig) {
			c.Illumina = types.IlluminaSpec{}
			c.PairedEnd = true
			c.TotalBatches = 9
		},
		"illumina binned paired region-free hmm": func(c *types.SequencerConfig) {
			c.Illumina = types.IlluminaSpec{}
			c.HMMModelPath = model
			c.PairedEnd = true
			c.TotalReads = 40
		},
		"illumina lognormal binned": func(c *types.SequencerConfig) {
			c.Illumina = types.IlluminaSpec{}
			c.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthLognormal, N50: 200, Median: 120}
			c.TotalBases = 9000
		},
		"illumina hmm paired lognormal": func(c *types.SequencerConfig) {
			c.HMMModelPath = model
			c.PairedEnd = true
			c.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthLognormal, N50: 200, Median: 120}
			c.TotalReads = 40
		},
		"ont gzip tier": func(c *types.SequencerConfig) {
			*c = scaleONTConfig(c.OutputDir)
			c.TotalBases = 60000
		},
		"ont hmm": func(c *types.SequencerConfig) {
			*c = scaleONTConfig(c.OutputDir)
			c.HMMModelPath = model
			c.TotalBases = 30000
		},
	}
	for name, edit := range cases {
		t.Run(name, func(t *testing.T) {
			c := goldenConfig(t.TempDir())
			c.Seed = 31
			edit(&c)
			c.Workers = 1
			one, e1 := runTree(t, c)
			c.OutputDir, c.Workers = t.TempDir(), 8
			eight, _ := runTree(t, c)
			if len(one) < 3 {
				t.Fatalf("only %d files written", len(one))
			}
			if _, ok := one[manifest.FileName]; !ok {
				t.Fatal("no manifest written")
			}
			if c.Platform == "ont" {
				if _, ok := one[nanopore.SummaryFileName]; !ok {
					t.Fatal("no sequencing summary written")
				}
			}
			for name, data := range one {
				if !bytes.Equal(data, eight[name]) {
					t.Errorf("%s differs between 1 and 8 workers", name)
				}
			}
			if len(one) != len(eight) {
				t.Errorf("%d files with 1 worker, %d with 8", len(one), len(eight))
			}
			if err := e1.Manifest().VerifyDir(e1.config.OutputDir); err != nil {
				t.Error(err)
			}
		})
	}
}

// fastqTree returns the parsed FASTQ records of every FASTQ file in tree,
// in file-name order.
func fastqTree(t *testing.T, tree map[string][]byte) map[string][]fastqRec {
	t.Helper()
	out := map[string][]fastqRec{}
	for name, data := range tree {
		if !strings.Contains(name, ".fastq") {
			continue
		}
		if strings.HasSuffix(name, ".gz") {
			p := filepath.Join(t.TempDir(), "x.gz")
			if err := os.WriteFile(p, data, 0o600); err != nil {
				t.Fatal(err)
			}
			data = readMaybeGzip(t, p)
		}
		out[name] = parseFASTQ(t, string(data))
	}
	return out
}

func sortedKeys[V any](m map[string]V) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// TotalBases lands within one read of the target on the actual (indel-
// changed) read lengths, files hold whole batches except the last, and the
// manifest totals agree with the files.
func TestScale_TotalBasesWithinOneRead(t *testing.T) {
	for _, workers := range []int{1, 3} {
		c := scaleONTConfig(t.TempDir())
		c.TotalBases = 100000
		c.Workers = workers
		tree, e := runTree(t, c)
		files := fastqTree(t, tree)
		var names []string // in output order
		for _, o := range e.Manifest().Outputs {
			if o.Role == manifest.RoleFASTQ {
				names = append(names, o.Path)
			}
		}
		if len(names) != len(files) || len(names) < 3 {
			t.Fatalf("manifest lists %d FASTQ files, output has %d", len(names), len(files))
		}
		var total, last int64
		for i, name := range names {
			recs := files[name]
			if i < len(names)-1 && len(recs) != c.ReadsPerBatch {
				t.Errorf("%s: %d reads, want a full batch of %d", name, len(recs), c.ReadsPerBatch)
			}
			for _, r := range recs {
				total += int64(len(r.seq))
				last = int64(len(r.seq))
			}
		}
		if total < c.TotalBases || total-last >= c.TotalBases {
			t.Errorf("workers %d: %d bases (last read %d) for a target of %d", workers, total, last, c.TotalBases)
		}
		m := e.Manifest()
		if m.Totals.Bases != total || m.FlowCells[0].Bases != total {
			t.Errorf("manifest bases %d / %d, files hold %d", m.Totals.Bases, m.FlowCells[0].Bases, total)
		}
		rows := readSummary(t, filepath.Join(c.OutputDir, nanopore.SummaryFileName))
		if int64(len(rows)) != m.Totals.Reads {
			t.Errorf("%d summary rows for %d reads", len(rows), m.Totals.Reads)
		}
	}

	// Paired-end: the target counts both mates; the last pair crosses it.
	c := goldenConfig(t.TempDir())
	c.Seed, c.PairedEnd, c.TotalBases = 4, true, 5000
	tree, _ := runTree(t, c)
	var total int64
	for _, recs := range fastqTree(t, tree) {
		for _, r := range recs {
			total += int64(len(r.seq))
		}
	}
	if total < 5000 || total-2*int64(c.ReadLength) >= 5000 {
		t.Errorf("paired: %d bases for a target of 5000", total)
	}
}

func TestScale_TotalReadsExact(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed, c.TotalReads, c.PairedEnd, c.Workers = 8, 20, true, 4
	tree, e := runTree(t, c)
	files := fastqTree(t, tree)
	want := map[string]int{
		"Gold_S1_L003_R1_001.fastq": 7, "Gold_S1_L003_R2_001.fastq": 7,
		"Gold_S1_L003_R1_002.fastq": 7, "Gold_S1_L003_R2_002.fastq": 7,
		"Gold_S1_L003_R1_003.fastq": 6, "Gold_S1_L003_R2_003.fastq": 6,
	}
	if len(files) != len(want) {
		t.Fatalf("files %v", sortedKeys(files))
	}
	for name, n := range want {
		if len(files[name]) != n {
			t.Errorf("%s: %d reads, want %d", name, len(files[name]), n)
		}
	}
	if m := e.Manifest(); m.Totals.Reads != 40 || m.Totals.Files != 6 {
		t.Errorf("manifest totals %+v", m.Totals)
	}
	// A partial batch is the prefix of the full batch.
	full := fileOutputs(t, func() types.SequencerConfig {
		d := goldenConfig("")
		d.Seed, d.PairedEnd, d.TotalBatches = 8, true, 3
		return d
	}())
	if !bytes.HasPrefix(full["Gold_S1_L003_R2_003.fastq"], tree["Gold_S1_L003_R2_003.fastq"]) {
		t.Error("the partial last file is not a prefix of the full batch")
	}

	// Streamed output honours the target too.
	c.TotalReads = 10
	if recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{}))); len(recs) != 20 {
		t.Errorf("stream: %d records, want 20 (10 pairs)", len(recs))
	}
}

func TestScale_StreamTotalBases(t *testing.T) {
	c := lognormalConfig("")
	c.TotalBatches, c.TotalBases = 0, 40000
	recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{})))
	ls := readLengths(recs)
	var total int
	for _, l := range ls {
		total += l
	}
	if total < 40000 || total-ls[len(ls)-1] >= 40000 {
		t.Errorf("stream: %d bases for a target of 40000", total)
	}
}

// Mates of a pair share their sampled length in file output, as in
// interleaved streams; the reads differ.
func TestScale_PairedFileMatesShareLength(t *testing.T) {
	c := lognormalConfig("")
	c.PairedEnd = true
	c.ReadsPerBatch = 200
	files := fileOutputs(t, c)
	for batch := 1; batch <= c.TotalBatches; batch++ {
		r1 := parseFASTQ(t, string(files["Gold_S1_L003_R1_00"+string(rune('0'+batch))+".fastq"]))
		r2 := parseFASTQ(t, string(files["Gold_S1_L003_R2_00"+string(rune('0'+batch))+".fastq"]))
		if len(r1) != c.ReadsPerBatch || len(r2) != len(r1) {
			t.Fatalf("batch %d: %d and %d reads", batch, len(r1), len(r2))
		}
		l1, l2 := readLengths(r1), readLengths(r2)
		if !slices.Equal(l1, l2) {
			t.Errorf("batch %d: mate lengths differ", batch)
		}
		if slices.Min(l1) == slices.Max(l1) {
			t.Errorf("batch %d: all reads have length %d", batch, l1[0])
		}
		if r1[0].seq == r2[0].seq {
			t.Errorf("batch %d: R1 and R2 reads are identical", batch)
		}
	}
}

// With an HMM model, a file generated alone equals that file of a full run:
// each file has its own model instance.
func TestScale_HMMFileIndependentOfOrder(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed, c.HMMModelPath, c.TotalBatches = 9, trainedModel(t), 3
	all := fileOutputs(t, c)
	e := NewEngine(c)
	defer e.Close()
	path := filepath.Join(c.OutputDir, "only3.fastq")
	if err := e.GenerateFASTQFile(path, 3, 1); err != nil {
		t.Fatal(err)
	}
	data, _ := os.ReadFile(path)
	if !bytes.Equal(data, all["Gold_S1_L003_003.fastq"]) {
		t.Error("batch 3 generated alone differs from batch 3 of a full run")
	}
}

func TestScale_FlowCells(t *testing.T) {
	c := scaleONTConfig(t.TempDir())
	c.FlowCells, c.TotalReads, c.Workers = 3, 30, 4
	tree, e := runTree(t, c)
	m := e.Manifest()
	if len(m.FlowCells) != 3 {
		t.Fatalf("manifest has %d flow cells", len(m.FlowCells))
	}
	runIDs, fcIDs, seeds := map[string]bool{}, map[string]bool{}, map[uint64]bool{}
	for _, fc := range m.FlowCells {
		runIDs[fc.RunID], fcIDs[fc.FlowCellID], seeds[fc.Seed] = true, true, true
		if fc.Dir != fc.FlowCellID || fc.Reads != 30 || fc.Files != 3 {
			t.Errorf("flow cell %+v: want dir = flow-cell ID, 30 reads in 3 files", fc)
		}
		var reads int
		for name, recs := range fastqTree(t, tree) {
			if !strings.HasPrefix(name, fc.Dir+"/") {
				continue
			}
			for _, r := range recs {
				m := ontHeader.FindStringSubmatch(r.header)
				if m == nil || m[2] != fc.RunID || m[6] != fc.FlowCellID {
					t.Fatalf("%s: header %q does not carry run %s / flow cell %s", name, r.header, fc.RunID, fc.FlowCellID)
				}
				reads++
			}
		}
		if reads != 30 {
			t.Errorf("flow cell %s: %d reads in its directory", fc.Dir, reads)
		}
		if rows := readSummary(t, filepath.Join(c.OutputDir, fc.Dir, nanopore.SummaryFileName)); len(rows) != 30 {
			t.Errorf("flow cell %s: %d summary rows", fc.Dir, len(rows))
		}
	}
	if len(runIDs) != 3 || len(fcIDs) != 3 || len(seeds) != 3 || seeds[e.Seed()] {
		t.Errorf("flow cells not independent: run IDs %v, flow cells %v, seeds %v", runIDs, fcIDs, seeds)
	}
	if m.Totals.Reads != 90 {
		t.Errorf("total reads %d", m.Totals.Reads)
	}

	// Illumina: explicit IDs get a per-flow-cell suffix.
	ic := goldenConfig(t.TempDir())
	ic.Seed, ic.FlowCells, ic.FlowCellID, ic.RunID = 2, 2, "HABCDEFXX", "run7"
	itree, ie := runTree(t, ic)
	for _, want := range []string{"HABCDEFXX-1/Gold_S1_L003_001.fastq", "HABCDEFXX-2/Gold_S1_L003_002.fastq", manifest.FileName} {
		if _, ok := itree[want]; !ok {
			t.Errorf("missing %s in %v", want, sortedKeys(itree))
		}
	}
	if fcs := ie.Manifest().FlowCells; fcs[1].RunID != "run7-2" || fcs[1].FlowCellID != "HABCDEFXX-2" {
		t.Errorf("flow cell 2 %+v", fcs[1])
	}
	if bytes.Equal(itree["HABCDEFXX-1/Gold_S1_L003_001.fastq"], itree["HABCDEFXX-2/Gold_S1_L003_001.fastq"]) {
		t.Error("flow cells have identical reads")
	}
}

func TestScale_MaxOutputCap(t *testing.T) {
	c := goldenConfig(t.TempDir())
	c.Seed, c.TotalBatches, c.MaxOutputGB = 5, 5, 1000e-9 // 1000 bases: 16 reads of 60
	tree, e := runTree(t, c)
	var total int
	for _, recs := range fastqTree(t, tree) {
		for _, r := range recs {
			total += len(r.seq)
		}
	}
	if total != 960 {
		t.Errorf("%d bases under a cap of 1000 (want 16 reads of 60)", total)
	}
	if fc := e.Manifest().FlowCells[0]; !fc.Capped || fc.Bases != 960 {
		t.Errorf("manifest flow cell %+v", fc)
	}
	if recs := parseFASTQ(t, string(streamOutput(t, c, GenerateOpts{}))); len(recs) != 16 {
		t.Errorf("stream: %d reads under the cap, want 16", len(recs))
	}

	c.OutputDir = t.TempDir()
	c.TotalBases = 2000
	en := NewEngine(c)
	defer en.Close()
	if err := en.Generate(); err == nil || !strings.Contains(err.Error(), "max output") {
		t.Errorf("total_bases above the cap: err = %v", err)
	}
}

func TestScale_Validation(t *testing.T) {
	cases := map[string]func(*types.SequencerConfig){
		"both targets":     func(c *types.SequencerConfig) { c.TotalBases, c.TotalReads = 100, 100 },
		"negative bases":   func(c *types.SequencerConfig) { c.TotalBases = -1 },
		"negative reads":   func(c *types.SequencerConfig) { c.TotalReads = -1 },
		"negative workers": func(c *types.SequencerConfig) { c.Workers = -2 },
		"ont shared flow cell id": func(c *types.SequencerConfig) {
			*c = scaleONTConfig(c.OutputDir)
			c.FlowCells, c.ONT.FlowCellID = 2, "PAW00001"
		},
		"bad s3 uri": func(c *types.SequencerConfig) { c.OutputDir = "s3://" },
		"pattern without file number": func(c *types.SequencerConfig) {
			c.TotalBatches, c.TotalBases, c.Output.FilePattern = 1, 1000, "{sample_name}.fastq"
		},
	}
	for name, edit := range cases {
		t.Run(name, func(t *testing.T) {
			c := goldenConfig(t.TempDir())
			edit(&c)
			e := NewEngine(c)
			defer e.Close()
			if err := e.Generate(); err == nil {
				t.Error("Generate accepted the configuration")
			}
			if err := e.GenerateToWriter(io.Discard, GenerateOpts{}); err == nil {
				t.Error("GenerateToWriter accepted the configuration")
			}
			if entries, _ := os.ReadDir(c.OutputDir); len(entries) != 0 {
				t.Errorf("output written despite the error: %d entries", len(entries))
			}
		})
	}
	c := goldenConfig(t.TempDir())
	c.FlowCells = 2
	e := NewEngine(c)
	defer e.Close()
	if err := e.GenerateToWriter(io.Discard, GenerateOpts{}); err == nil || !strings.Contains(err.Error(), "file output") {
		t.Errorf("flow cells on a stream: err = %v", err)
	}
	// One file when the reads target fits in a batch: no {file_number} needed.
	c = goldenConfig(t.TempDir())
	c.TotalBatches, c.TotalReads, c.Output.FilePattern = 3, 5, "{sample_name}.fastq"
	tree, _ := runTree(t, c)
	if len(parseFASTQ(t, string(tree["Gold.fastq"]))) != 5 {
		t.Errorf("single-file run: %v", sortedKeys(tree))
	}
}

func TestManifest_MatchesFilesAndOmitsKeys(t *testing.T) {
	dir := t.TempDir()
	keyPath := filepath.Join(t.TempDir(), "secret.key")
	const key = "correct-horse-battery-staple-0123456789" // gitleaks:allow -- test fixture, not a real key
	if err := os.WriteFile(keyPath, []byte(key+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	c := scaleONTConfig(dir)
	c.Fingerprinting = types.FingerprintConfig{Enable: true, CustomerFingerprint: true, SteganographyLayers: []string{"read_ids"}, KeyFile: keyPath}
	c.Enterprise.CustomerID = "cust-9"
	c.Enterprise.Fingerprinting.KeyFile = keyPath
	tree, e := runTree(t, c)
	raw := tree[manifest.FileName]
	for _, secret := range []string{key, keyPath, "secret.key", dir} {
		if bytes.Contains(raw, []byte(secret)) {
			t.Errorf("manifest contains %q", secret)
		}
	}
	m, err := manifest.Read(filepath.Join(dir, manifest.FileName))
	if err != nil {
		t.Fatal(err)
	}
	if err := m.VerifyDir(dir); err != nil {
		t.Error(err)
	}
	if m.Seed != e.Seed() || m.Config.Seed != e.Seed() || m.Tool.Name != "phantom-science" || m.Config.Workers != 0 || m.Config.OutputDir != "" {
		t.Errorf("manifest header %+v", m)
	}
	if !m.RunStart.Equal(e.ontRunInfo().Start) {
		t.Errorf("run_start %v, want the seed-derived ONT start %v", m.RunStart, e.ontRunInfo().Start)
	}
	roles := map[string]int{}
	var reads int64
	for _, o := range m.Outputs {
		roles[o.Role]++
		reads += o.Reads
		if _, ok := tree[o.Path]; !ok {
			t.Errorf("manifest lists %s, not in the output", o.Path)
		}
	}
	if roles[manifest.RoleFASTQ] != len(e.GeneratedFiles()) || roles[manifest.RoleSummary] != 1 || len(tree) != len(m.Outputs)+1 {
		t.Errorf("roles %v for %d files in the output", roles, len(tree))
	}
	if reads != m.Totals.Reads || reads != int64(c.TotalBatches*c.ReadsPerBatch) {
		t.Errorf("manifest reads %d, totals %d", reads, m.Totals.Reads)
	}
	var generic map[string]any
	if err := json.Unmarshal(raw, &generic); err != nil {
		t.Fatal(err)
	}
	for _, k := range []string{"format", "version", "tool", "seed", "run_start", "config", "flow_cells", "outputs", "totals"} {
		if _, ok := generic[k]; !ok {
			t.Errorf("manifest lacks %q", k)
		}
	}
}

func sha256Hex(b []byte) string {
	h := sha256.Sum256(b)
	return hex.EncodeToString(h[:])
}

// fakeUploader records uploads in memory and checks the spool file.
type fakeUploader struct {
	mu      sync.Mutex
	objects map[string][]byte
	tags    map[string]map[string]string
	spools  []string
}

func (f *fakeUploader) Upload(_ context.Context, in s3sink.UploadInput) error {
	data, err := io.ReadAll(in.Body)
	if err != nil {
		return err
	}
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.objects == nil {
		f.objects, f.tags = map[string][]byte{}, map[string]map[string]string{}
	}
	key := in.Bucket + "/" + in.Key
	if _, dup := f.objects[key]; dup {
		return errors.New("object uploaded twice: " + key)
	}
	f.objects[key], f.tags[key] = data, in.Tags
	f.spools = append(f.spools, in.Body.Name())
	return nil
}

// s3Run generates c to a fake S3 prefix and returns the uploads.
func s3Run(t *testing.T, c types.SequencerConfig) (*fakeUploader, *Engine) {
	t.Helper()
	up := &fakeUploader{}
	e := NewEngine(c)
	t.Cleanup(e.Close)
	e.SetUploader(up)
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	for _, s := range up.spools {
		if _, err := os.Stat(filepath.Dir(s)); !os.IsNotExist(err) {
			t.Errorf("spool directory %s not removed", filepath.Dir(s))
		}
	}
	return up, e
}

func TestS3_UploadsMatchManifestAndLocal(t *testing.T) {
	c := scaleONTConfig(t.TempDir())
	c.FlowCells, c.TotalBases, c.Workers = 2, 40000, 4
	local, _ := runTree(t, c)

	c.OutputDir = "s3://bucket/runs/r1"
	up, e := s3Run(t, c)
	if len(up.objects) != len(local) {
		t.Errorf("%d objects uploaded, %d local files", len(up.objects), len(local))
	}
	for rel, data := range local {
		got, ok := up.objects["bucket/runs/r1/"+rel]
		if !ok {
			t.Errorf("no object for %s", rel)
			continue
		}
		if !bytes.Equal(got, data) && rel != manifest.FileName {
			t.Errorf("%s: uploaded bytes differ from local output", rel)
		}
	}
	m := e.Manifest()
	var mf manifest.Manifest
	if err := json.Unmarshal(up.objects["bucket/runs/r1/"+manifest.FileName], &mf); err != nil {
		t.Fatal(err)
	}
	if mf.Destination != "s3://bucket/runs/r1/" || m.Destination != mf.Destination {
		t.Errorf("destination %q", mf.Destination)
	}
	var lm manifest.Manifest
	if err := json.Unmarshal(local[manifest.FileName], &lm); err != nil {
		t.Fatal(err)
	}
	for i, o := range mf.Outputs {
		data, ok := up.objects["bucket/"+o.Key]
		if !ok || o.Key != "runs/r1/"+o.Path {
			t.Errorf("output %s: key %q not uploaded", o.Path, o.Key)
			continue
		}
		if int64(len(data)) != o.Bytes || sha256Hex(data) != o.SHA256 {
			t.Errorf("%s: manifest size/sha256 do not match the uploaded object", o.Key)
		}
		if lo := lm.Outputs[i]; lo.Path != o.Path || lo.SHA256 != o.SHA256 {
			t.Errorf("output %d differs between the local and S3 manifests", i)
		}
		if tags := up.tags["bucket/"+o.Key]; tags["phantom:regenerable"] != "seed" || tags["phantom:output-bytes"] == "" {
			t.Errorf("%s: tags %v", o.Key, tags)
		}
	}
	if gf := e.GeneratedFiles(); len(gf) == 0 || !strings.HasPrefix(gf[0], "s3://bucket/runs/r1/") {
		t.Errorf("GeneratedFiles %v", gf)
	}
}
