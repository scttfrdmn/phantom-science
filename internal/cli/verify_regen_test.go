//go:build !release

package cli

import (
	"bytes"
	"compress/gzip"
	"encoding/json"
	"io"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Verification of keyed-seed fingerprints by regeneration: suspect reads are
// matched against the candidate runs regenerated from the key and the run's
// records. These tests need no stegano library.

const regenKey = "regen-verify-key-0123456789abcdef" // gitleaks:allow -- test fixture

func regenKeyFile(t *testing.T, key string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "k.key")
	if err := os.WriteFile(p, []byte(key+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

// illuminaRun returns a binned paired-end config: three batches of 200 pairs.
func illuminaRun(dir string) types.SequencerConfig {
	c := testSequencerConfig(dir)
	c.PairedEnd = true
	c.TotalBatches = 3
	c.Output.Compression = "gzip"
	return c
}

func ontRun(dir string) types.SequencerConfig {
	c := testSequencerConfig(dir)
	c.Platform = "ont"
	c.ReadLength = 2000
	c.ReadsPerBatch = 60
	c.TotalBatches = 2
	c.ErrorModel.Tier = "hac"
	return c
}

// keyedRun generates c fingerprinted for customer/dataset with salt and
// returns its files and manifest path.
func keyedRun(t *testing.T, c types.SequencerConfig, keyFile, customer, dataset string, salt uint64) ([]string, string) {
	t.Helper()
	c.Enterprise.CustomerID, c.DatasetID, c.Seed = customer, dataset, salt
	c.Fingerprinting = types.FingerprintConfig{Enable: true, KeyFile: keyFile}
	return generateFASTQ(t, c), filepath.Join(c.OutputDir, manifest.FileName)
}

type rec struct{ head, seq, qual string }

func readRecs(t *testing.T, path string) []rec {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		zr, err := gzip.NewReader(f)
		if err != nil {
			t.Fatal(err)
		}
		r = zr
	}
	data, err := io.ReadAll(r)
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSuffix(string(data), "\n"), "\n")
	var out []rec
	for i := 0; i+3 < len(lines); i += 4 {
		out = append(out, rec{lines[i], lines[i+1], lines[i+3]})
	}
	return out
}

func writeRecs(t *testing.T, recs []rec) string {
	t.Helper()
	var b strings.Builder
	for _, r := range recs {
		b.WriteString(r.head + "\n" + r.seq + "\n+\n" + r.qual + "\n")
	}
	p := filepath.Join(t.TempDir(), "suspect.fastq")
	if err := os.WriteFile(p, []byte(b.String()), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

func verifyOut(t *testing.T, paths []string, opts VerifyOptions) (string, error) {
	t.Helper()
	var out bytes.Buffer
	err := VerifyFiles(&out, paths, opts)
	return out.String(), err
}

func TestVerifyRegeneration(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := regenKeyFile(t, regenKey)
	filesA, manA := keyedRun(t, illuminaRun(t.TempDir()), keyFile, "ACME", "study_1", 11)
	_, manB := keyedRun(t, illuminaRun(t.TempDir()), keyFile, "BETA", "study_1", 11)
	both := VerifyOptions{KeyFile: keyFile, Manifests: []string{manB, manA}}

	third := readRecs(t, filesA[2]) // batch 2, R1
	rng := rand.New(rand.NewPCG(1, 2))
	shuffled := append([]rec(nil), third...)
	rng.Shuffle(len(shuffled), func(i, j int) { shuffled[i], shuffled[j] = shuffled[j], shuffled[i] })
	var requal, trimmed []rec
	for _, r := range third {
		requal = append(requal, rec{"@x", r.seq, strings.Repeat("F", len(r.qual))})
		trimmed = append(trimmed, rec{r.head, r.seq[10 : len(r.seq)-10], r.qual[10 : len(r.qual)-10]})
	}
	random := make([]rec, 50)
	for i := range random {
		s := make([]byte, 150)
		for j := range s {
			s[j] = "ACGT"[rng.IntN(4)]
		}
		random[i] = rec{"@r", string(s), strings.Repeat("F", 150)}
	}
	unrelated := generateFASTQ(t, illuminaRun(t.TempDir()))[0] // unfingerprinted, seed drawn

	positive := map[string]string{
		"whole file":           filesA[0],
		"20 reads of file 3":   writeRecs(t, third[50:70]),
		"shuffled":             writeRecs(t, shuffled),
		"qualities all F":      writeRecs(t, requal),
		"trimmed 10 each side": writeRecs(t, trimmed),
	}
	for name, path := range positive {
		t.Run(name, func(t *testing.T) {
			out, err := verifyOut(t, []string{path}, both)
			if err != nil {
				t.Fatalf("%v\n%s", err, out)
			}
			if !strings.Contains(out, "✅") || !strings.Contains(out, "matches "+filepath.Base(filepath.Dir(manA))) || !strings.Contains(out, "customer=ACME") {
				t.Errorf("not attributed to ACME:\n%s", out)
			}
			if strings.Contains(out, "also matches") {
				t.Errorf("matched a second candidate:\n%s", out)
			}
		})
	}
	// Trimmed reads match by containment, not exactly.
	if out, _ := verifyOut(t, []string{positive["trimmed 10 each side"]}, VerifyOptions{KeyFile: keyFile, Manifests: []string{manA}}); !strings.Contains(out, "200 of 200 reads (0 exact, 200 by 31-mer containment)") {
		t.Errorf("trimmed reads:\n%s", out)
	}

	negative := map[string]struct {
		paths []string
		opts  VerifyOptions
	}{
		"customer B only": {[]string{filesA[0]}, VerifyOptions{KeyFile: keyFile, Manifests: []string{manB}}},
		"wrong key":       {[]string{filesA[0]}, VerifyOptions{KeyFile: regenKeyFile(t, "another-key-0123456789abcdef"), Manifests: []string{manA}}},
		"unrelated run":   {[]string{unrelated}, both},
		"random reads":    {[]string{writeRecs(t, random)}, both},
	}
	for name, tc := range negative {
		t.Run(name, func(t *testing.T) {
			out, err := verifyOut(t, tc.paths, tc.opts)
			if err == nil || strings.Contains(out, "✅") {
				t.Fatalf("matched: %v\n%s", err, out)
			}
			if !strings.Contains(out, "❌") || !strings.Contains(out, " 0 of ") {
				t.Errorf("want a no-match verdict with 0 matched reads:\n%s", out)
			}
		})
	}

	// One verdict per file.
	out, err := verifyOut(t, []string{filesA[1], unrelated}, both)
	if err == nil || strings.Count(out, "✅") != 1 || strings.Count(out, "❌") != 1 {
		t.Errorf("mixed files: %v\n%s", err, out)
	}

	// Limits: regenerating only the first file cannot attribute reads of
	// the third; MaxReads caps the suspect reads examined.
	if out, err := verifyOut(t, []string{positive["20 reads of file 3"]}, VerifyOptions{KeyFile: keyFile, Manifests: []string{manA}, MaxFiles: 1}); err == nil {
		t.Errorf("MaxFiles 1 matched reads of file 3:\n%s", out)
	}
	if out, _ := verifyOut(t, []string{filesA[0]}, VerifyOptions{KeyFile: keyFile, Manifests: []string{manA}, MaxReads: 25}); !strings.Contains(out, "25 of 25 reads") {
		t.Errorf("MaxReads 25:\n%s", out)
	}

	// A decoy run is measured and matches nothing.
	if out, err := verifyOut(t, []string{filesA[0]}, VerifyOptions{KeyFile: keyFile, Manifests: []string{manA}, Decoy: true}); err != nil || !strings.Contains(out, "decoy run 0") {
		t.Errorf("decoy: %v\n%s", err, out)
	}
}

// --config with --customer/--dataset/--salt, and a candidates file, select
// the run without a manifest.
func TestVerifyRegeneration_ConfigAndCandidates(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := regenKeyFile(t, regenKey)
	dir := t.TempDir()
	c := ontRun(filepath.Join(dir, "out"))
	files, _ := keyedRun(t, c, keyFile, "ACME", "study_9", 7)

	cfgPath := filepath.Join(dir, "run.json")
	data, err := json.Marshal(c)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(cfgPath, data, 0o600); err != nil {
		t.Fatal(err)
	}
	out, err := verifyOut(t, files, VerifyOptions{KeyFile: keyFile, ConfigPath: cfgPath, Customer: "ACME", Dataset: "study_9", Salt: 7, SaltSet: true})
	if err != nil || strings.Count(out, "✅") != len(files) {
		t.Fatalf("--config: %v\n%s", err, out)
	}
	if out, err := verifyOut(t, files[:1], VerifyOptions{KeyFile: keyFile, ConfigPath: cfgPath, Customer: "ACME", Dataset: "study_9", Salt: 8, SaltSet: true}); err == nil {
		t.Errorf("wrong salt matched:\n%s", out)
	}

	cands := filepath.Join(dir, "candidates.tsv")
	tsv := "customer\tdataset\tsalt\n# suspects\nBETA\tstudy_9\t7\nACME\tstudy_9\t7\nACME\tstudy_8\t7\n"
	if err := os.WriteFile(cands, []byte(tsv), 0o600); err != nil {
		t.Fatal(err)
	}
	trimmed := []rec{}
	for _, r := range readRecs(t, files[1]) {
		trimmed = append(trimmed, rec{"@read", r.seq[10 : len(r.seq)-10], r.qual[10 : len(r.qual)-10]})
	}
	out, err = verifyOut(t, []string{writeRecs(t, trimmed)}, VerifyOptions{KeyFile: keyFile, ConfigPath: cfgPath, Candidates: cands})
	if err != nil || !strings.Contains(out, "✅") || !strings.Contains(out, "matches customer=ACME dataset=study_9 salt=7") {
		t.Fatalf("candidates: %v\n%s", err, out)
	}
	if strings.Count(out, "regenerated") != 3 {
		t.Errorf("want three candidates regenerated:\n%s", out)
	}
}

func TestVerifyRegeneration_Errors(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := regenKeyFile(t, regenKey)
	files, man := keyedRun(t, illuminaRun(t.TempDir()), keyFile, "ACME", "study_1", 11)
	for name, tc := range map[string]struct {
		opts VerifyOptions
		want string
	}{
		"missing manifest":          {VerifyOptions{KeyFile: keyFile, Manifests: []string{filepath.Join(t.TempDir(), manifest.FileName)}}, "reading manifest"},
		"no key":                    {VerifyOptions{Manifests: []string{man}}, "needs the fingerprint key"},
		"config without candidate":  {VerifyOptions{KeyFile: keyFile, ConfigPath: man}, "--customer"},
		"manifest and customer":     {VerifyOptions{KeyFile: keyFile, Manifests: []string{man}, Customer: "ACME"}, "--manifest records"},
		"candidates without config": {VerifyOptions{KeyFile: keyFile, Candidates: man}, "needs the run configuration"},
	} {
		t.Run(name, func(t *testing.T) {
			_, err := verifyOut(t, files[:1], tc.opts)
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Errorf("err = %v, want %q", err, tc.want)
			}
		})
	}
}

// An unfingerprinted run's manifest is a candidate too (by its seed), and
// needs no key.
func TestVerifyRegeneration_UnfingerprintedManifest(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	c := illuminaRun(t.TempDir())
	c.Seed = 3
	files := generateFASTQ(t, c)
	out, err := verifyOut(t, files[:2], VerifyOptions{Manifests: []string{filepath.Join(c.OutputDir, manifest.FileName)}})
	if err != nil || strings.Count(out, "✅") != 2 || !strings.Contains(out, "seed=3 (not fingerprinted)") {
		t.Errorf("%v\n%s", err, out)
	}
}

// A keyed run generated with --variants (#32) is attributed when verify is
// given the same VCF: every read is reproduced exactly. Without the VCF, or
// with another one, verify fails naming the VCF the manifest records,
// rather than regenerating the run without its variants.
func TestVerifyRegeneration_Variants(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := regenKeyFile(t, regenKey)
	dir := t.TempDir()
	writeVCF := func(name, body string) string {
		p := filepath.Join(dir, name)
		if err := os.WriteFile(p, []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+body), 0o600); err != nil {
			t.Fatal(err)
		}
		return p
	}
	vcf := writeVCF("calls.vcf", "chr5\t1000\t.\tA\tAGGG\t.\t.\tAF=0.5\nchr5\t1100\t.\tC\tT\t.\t.\tAF=0.5\n")
	other := writeVCF("other.vcf", "chr5\t1000\t.\tA\tAGGG\t.\t.\tAF=0.4\nchr5\t1100\t.\tC\tT\t.\t.\tAF=0.5\n")
	c := illuminaRun(t.TempDir())
	c.VariantsPath = vcf
	c.VariantsRegion = "chr5:1-3000" // wide enough for the paired fragments
	files, man := keyedRun(t, c, keyFile, "ACME", "study_v", 5)

	out, err := verifyOut(t, files[:2], VerifyOptions{KeyFile: keyFile, Manifests: []string{man}, Variants: []string{other, vcf}})
	if err != nil || strings.Count(out, "✅") != 2 || !strings.Contains(out, "200 of 200 reads (200 exact") {
		t.Errorf("with the VCF: %v\n%s", err, out)
	}
	for name, vs := range map[string][]string{"no VCF": nil, "another VCF": {other}} {
		out, err := verifyOut(t, files[:1], VerifyOptions{KeyFile: keyFile, Manifests: []string{man}, Variants: vs})
		if err == nil || !strings.Contains(err.Error(), "calls.vcf") || !strings.Contains(err.Error(), "--variants") {
			t.Errorf("%s: err = %v\n%s", name, err, out)
		}
	}
}
