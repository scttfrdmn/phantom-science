//go:build !release

package genomics

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Keyed-seed fingerprints (docs/CHARTER.md, principle 2): a fingerprinted
// run changes nothing but the seed.

const keyedTestKey = "keyed-seed-test-key-0123456789abcdef" // gitleaks:allow -- test fixture

// writeKeyFile writes key to a mode 0600 file and returns its path.
func writeKeyFile(t *testing.T, key string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "fingerprint.key")
	if err := os.WriteFile(p, []byte(key+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

// fingerprinted returns c with a keyed-seed fingerprint for customer and
// dataset, keyed by keyFile.
func fingerprinted(c types.SequencerConfig, keyFile, customer, dataset string) types.SequencerConfig {
	c.Enterprise.CustomerID = customer
	c.DatasetID = dataset
	c.Fingerprinting = types.FingerprintConfig{Enable: true, KeyFile: keyFile}
	return c
}

func deriveSeed(t *testing.T, key, customer, dataset string, salt uint64) uint64 {
	t.Helper()
	s, err := fingerprinting.DeriveSeed([]byte(key), customer, dataset, salt)
	if err != nil {
		t.Fatal(err)
	}
	return s
}

// withoutManifest returns tree's output files, without the manifest.
func withoutManifest(tree map[string][]byte) map[string][]byte {
	out := map[string][]byte{}
	for k, v := range tree {
		if k != manifest.FileName {
			out[k] = v
		}
	}
	return out
}

func sameTree(t *testing.T, what string, got, want map[string][]byte) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s: %d files %v, want %d %v", what, len(got), sortedKeys(got), len(want), sortedKeys(want))
	}
	for name, w := range want {
		if g, ok := got[name]; !ok {
			t.Errorf("%s: %s missing", what, name)
		} else if !bytes.Equal(g, w) {
			t.Errorf("%s: %s differs (%d vs %d bytes)", what, name, len(g), len(w))
		}
	}
}

// A fingerprinted run is byte-identical, file for file, to an
// unfingerprinted run whose seed is the derived seed, for binned paired
// Illumina (gzip, several batches) and ONT (two flow cells), with one or
// several workers.
func TestKeyedSeed_OutputEqualsUnfingerprintedRun(t *testing.T) {
	keyFile := writeKeyFile(t, keyedTestKey)
	illumina := func(dir string) types.SequencerConfig {
		c := binnedConfig(dir)
		c.ReadsPerBatch, c.TotalBatches = 300, 3
		c.Output.Compression = "gzip"
		return c
	}
	ont := func(dir string) types.SequencerConfig {
		c := scaleONTConfig(dir)
		c.FlowCells = 2
		return c
	}
	for name, mk := range map[string]func(string) types.SequencerConfig{"illumina_paired_binned": illumina, "ont": ont} {
		t.Run(name, func(t *testing.T) {
			const salt = 4242
			derived := deriveSeed(t, keyedTestKey, "ACME", "study_1", salt)

			plainCfg := mk(t.TempDir())
			plainCfg.Seed = derived
			plain, _ := runTree(t, plainCfg)

			for _, workers := range []int{1, 3} {
				c := fingerprinted(mk(t.TempDir()), keyFile, "ACME", "study_1")
				c.Seed, c.Workers = salt, workers
				fp, e := runTree(t, c)
				if e.Seed() != derived {
					t.Fatalf("engine seed %d, want the derived seed %d", e.Seed(), derived)
				}
				sameTree(t, "workers="+strconv.Itoa(workers), withoutManifest(fp), withoutManifest(plain))
			}
		})
	}
}

// qualityAlphabet returns the set of quality characters in tree's FASTQ.
func qualityAlphabet(t *testing.T, tree map[string][]byte) string {
	t.Helper()
	seen := map[rune]bool{}
	for _, recs := range fastqTree(t, tree) {
		for _, r := range recs {
			for _, q := range r.qual {
				seen[q] = true
			}
		}
	}
	var b strings.Builder
	for c := rune('!'); c <= '~'; c++ {
		if seen[c] {
			b.WriteRune(c)
		}
	}
	return b.String()
}

// Fingerprinted binned output uses exactly the four RTA3 bins, like any
// unfingerprinted run (the quality-score layer moved values off-bin).
func TestKeyedSeed_BinnedQualitiesStayOnBins(t *testing.T) {
	keyFile := writeKeyFile(t, keyedTestKey)
	c := fingerprinted(binnedConfig(t.TempDir()), keyFile, "ACME", "study_1")
	c.ReadsPerBatch = 2000
	fp, _ := runTree(t, c)
	pc := binnedConfig(t.TempDir())
	pc.ReadsPerBatch = c.ReadsPerBatch
	plain, _ := runTree(t, pc)
	got, want := qualityAlphabet(t, fp), qualityAlphabet(t, plain)
	// Q2 no-calls are rare, so whether a run of this size has a '#' is
	// chance (the two runs have different seeds); compare the other bins.
	noQ2 := func(s string) string { return strings.ReplaceAll(s, "#", "") }
	if noQ2(got) != noQ2(want) || strings.Trim(got, binChars) != "" || !strings.Contains(got, "F") {
		t.Errorf("quality alphabet %q, unfingerprinted %q, want the same subset of the bins %q", got, want, binChars)
	}
}

// The manifest records the scheme, customer, dataset and salt, and the
// configured seed is the salt; it holds neither the derived seed (of the run
// or of a flow cell) nor the key.
func TestKeyedSeed_ManifestRecordsSaltNotSeed(t *testing.T) {
	keyFile := writeKeyFile(t, keyedTestKey)
	dir := t.TempDir()
	c := fingerprinted(scaleONTConfig(dir), keyFile, "ACME", "study_1")
	c.Seed, c.FlowCells = 99, 2
	tree, e := runTree(t, c)
	raw := tree[manifest.FileName]
	derived := e.Seed()
	secrets := []string{keyedTestKey, keyFile, strconv.FormatUint(derived, 10)}
	for k := 1; k <= 2; k++ {
		for attempt := 0; attempt < 4; attempt++ {
			secrets = append(secrets, strconv.FormatUint(flowCellSeed(derived, k, attempt), 10))
		}
	}
	for _, s := range secrets {
		if bytes.Contains(raw, []byte(s)) {
			t.Errorf("manifest contains %q", s)
		}
	}
	m, err := manifest.Read(filepath.Join(dir, manifest.FileName))
	if err != nil {
		t.Fatal(err)
	}
	want := manifest.Fingerprint{Scheme: fingerprinting.KeyedSeedScheme, Customer: "ACME", Dataset: "study_1", Salt: 99}
	if m.Fingerprint == nil || *m.Fingerprint != want {
		t.Errorf("fingerprint %+v, want %+v", m.Fingerprint, want)
	}
	if m.Seed != 0 || m.Config.Seed != 99 {
		t.Errorf("seed %d, config seed %d; want 0 and the salt 99", m.Seed, m.Config.Seed)
	}
	for _, fc := range m.FlowCells {
		if fc.Seed != 0 {
			t.Errorf("flow cell %d records seed %d", fc.Index, fc.Seed)
		}
	}
	if !bytes.Contains(raw, []byte(`"scheme": "keyed_seed/v1"`)) || !bytes.Contains(raw, []byte(`"salt": "99"`)) {
		t.Errorf("manifest JSON lacks the scheme or salt:\n%s", raw)
	}
	if err := m.VerifyDir(dir); err != nil {
		t.Error(err)
	}
}

// Without --seed the salt is drawn, logged and recorded; generating again
// with it as the seed and the same key reproduces the run.
func TestKeyedSeed_DrawnSaltReproduces(t *testing.T) {
	keyFile := writeKeyFile(t, keyedTestKey)
	c := fingerprinted(binnedConfig(t.TempDir()), keyFile, "ACME", "study_1")
	c.Seed = 0
	var log bytes.Buffer
	e := NewEngineLog(c, &log)
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	salt, ok := e.Fingerprinted()
	if !ok || salt == 0 {
		t.Fatalf("Fingerprinted() = %d, %v", salt, ok)
	}
	if !strings.Contains(log.String(), "fingerprint salt "+strconv.FormatUint(salt, 10)) {
		t.Errorf("salt %d not logged:\n%s", salt, log.String())
	}
	if strings.Contains(log.String(), strconv.FormatUint(e.Seed(), 10)) {
		t.Errorf("derived seed logged:\n%s", log.String())
	}
	m, err := manifest.Read(filepath.Join(c.OutputDir, manifest.FileName))
	if err != nil {
		t.Fatal(err)
	}
	if m.Fingerprint == nil || m.Fingerprint.Salt != salt || m.Config.Seed != salt {
		t.Fatalf("manifest fingerprint %+v config seed %d, want salt %d", m.Fingerprint, m.Config.Seed, salt)
	}
	first, _ := runTree(t, func() types.SequencerConfig { c2 := c; c2.OutputDir = t.TempDir(); c2.Seed = salt; return c2 }())
	again := map[string][]byte{}
	for _, f := range e.GeneratedFiles() {
		data, err := os.ReadFile(f)
		if err != nil {
			t.Fatal(err)
		}
		again[filepath.Base(f)] = data
	}
	sameTree(t, "reproduced from the salt", withoutManifest(first), again)
}

// Streamed output with a keyed seed equals streamed output with the derived
// seed.
func TestKeyedSeed_StreamEqualsDerivedSeed(t *testing.T) {
	keyFile := writeKeyFile(t, keyedTestKey)
	c := fingerprinted(binnedConfig(""), keyFile, "ACME", "study_1")
	c.Seed = 5
	var fp, plain bytes.Buffer
	if err := NewEngineLog(c, nil).GenerateToWriter(&fp, GenerateOpts{}); err != nil {
		t.Fatal(err)
	}
	p := binnedConfig("")
	p.Seed = deriveSeed(t, keyedTestKey, "ACME", "study_1", 5)
	if err := NewEngineLog(p, nil).GenerateToWriter(&plain, GenerateOpts{}); err != nil {
		t.Fatal(err)
	}
	if fp.Len() == 0 || !bytes.Equal(fp.Bytes(), plain.Bytes()) {
		t.Errorf("streamed keyed-seed output (%d bytes) differs from the derived-seed run (%d bytes)", fp.Len(), plain.Len())
	}
}

// Unusable fingerprint settings fail before anything is written.
func TestKeyedSeed_ConfigErrors(t *testing.T) {
	keyFile := writeKeyFile(t, keyedTestKey)
	t.Setenv(fingerprinting.KeyEnvVar, "")
	cases := map[string]struct {
		edit func(*types.SequencerConfig)
		want string
	}{
		"no key":         {func(c *types.SequencerConfig) { c.Fingerprinting.KeyFile = "" }, "needs a key"},
		"no customer":    {func(c *types.SequencerConfig) { c.Enterprise.CustomerID = "" }, "customer"},
		"no dataset":     {func(c *types.SequencerConfig) { c.DatasetID = "" }, "dataset"},
		"quality_scores": {func(c *types.SequencerConfig) { c.Fingerprinting.SteganographyLayers = []string{"quality_scores"} }, "CHARTER"},
		"unknown layer":  {func(c *types.SequencerConfig) { c.Fingerprinting.SteganographyLayers = []string{"headers"} }, "unknown"},
		"exposed key": {func(c *types.SequencerConfig) {
			p := writeKeyFile(t, keyedTestKey)
			os.Chmod(p, 0o644) // #nosec G302 -- the test needs an exposed key file
			c.Fingerprinting.KeyFile = p
		}, "chmod 600"},
	}
	for name, tc := range cases {
		t.Run(name, func(t *testing.T) {
			dir := t.TempDir()
			c := fingerprinted(binnedConfig(dir), keyFile, "ACME", "study_1")
			tc.edit(&c)
			e := NewEngineLog(c, nil)
			err := e.Generate()
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("Generate: %v, want an error containing %q", err, tc.want)
			}
			if name == "quality_scores" && !errors.Is(err, fingerprinting.ErrQualityScoresLayer) {
				t.Errorf("error %v is not ErrQualityScoresLayer", err)
			}
			entries, _ := os.ReadDir(dir)
			if len(entries) != 0 {
				t.Errorf("files written before the error: %v", entries)
			}
			if err := e.GenerateToWriter(&bytes.Buffer{}, GenerateOpts{}); err == nil {
				t.Error("GenerateToWriter accepted the configuration")
			}
		})
	}
}

// The key may come from the environment.
func TestKeyedSeed_KeyFromEnvironment(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, keyedTestKey)
	c := fingerprinted(binnedConfig(t.TempDir()), "", "ACME", "study_1")
	e := NewEngineLog(c, nil)
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	if want := deriveSeed(t, keyedTestKey, "ACME", "study_1", c.Seed); e.Seed() != want {
		t.Errorf("seed %d, want %d derived with the environment key", e.Seed(), want)
	}
}

// The read_ids layer alone writes the labels and leaves the seed alone.
func TestKeyedSeed_ReadIDsLayerAloneHasNoKeyedSeed(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	c := binnedConfig(t.TempDir())
	c.Enterprise.CustomerID = "ACME"
	c.Fingerprinting = types.FingerprintConfig{Enable: true, CustomerFingerprint: true, SteganographyLayers: []string{"read_ids"}}
	e := NewEngineLog(c, nil)
	if err := e.Generate(); err != nil {
		t.Fatal(err)
	}
	if _, ok := e.Fingerprinted(); ok || e.Seed() != c.Seed {
		t.Errorf("read_ids layer derived a seed (%d)", e.Seed())
	}
}

// Uploaded objects of a keyed-seed run are tagged as regenerable with the
// key and carry no seed.
func TestKeyedSeed_S3TagsCarryNoSeed(t *testing.T) {
	keyFile := writeKeyFile(t, keyedTestKey)
	c := fingerprinted(scaleONTConfig(""), keyFile, "ACME", "study_1")
	c.OutputDir = "s3://bucket/runs/fp"
	up, e := s3Run(t, c)
	if len(up.tags) == 0 {
		t.Fatal("nothing uploaded")
	}
	for key, tags := range up.tags {
		if strings.HasSuffix(key, manifest.FileName) {
			continue // the manifest has no tags
		}
		if _, ok := tags["phantom:seed"]; ok || tags["phantom:regenerable"] != "keyed_seed" {
			t.Errorf("%s: tags %v", key, tags)
		}
		for _, v := range tags {
			if v == strconv.FormatUint(e.Seed(), 10) {
				t.Errorf("%s: a tag holds the derived seed", key)
			}
		}
	}
}
