//go:build !release

package cli

import (
	"bytes"
	"encoding/json"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Input files of verified runs (#37): a manifest records the reference
// FASTA, HMM model and read-length histogram by name, size and sha256, not
// by path, and verify regenerates the run only with the files that match.

// randomBases returns n uniform bases from seed.
func randomBases(seed uint64, n int) []byte {
	r := rand.New(rand.NewPCG(seed, seed))
	b := make([]byte, n)
	for i := range b {
		b[i] = "ACGT"[r.IntN(4)]
	}
	return b
}

// writeFASTA writes contigs chr1 and chr2 of 20 kb drawn from seed into dir.
func writeFASTA(t *testing.T, dir, name string, seed uint64) (string, map[string][]byte) {
	t.Helper()
	seqs := map[string][]byte{"chr1": randomBases(seed, 20000), "chr2": randomBases(seed+1, 20000)}
	var buf bytes.Buffer
	for _, c := range []string{"chr1", "chr2"} {
		buf.WriteString(">" + c + "\n")
		for i := 0; i < len(seqs[c]); i += 60 {
			buf.Write(seqs[c][i:min(i+60, len(seqs[c]))])
			buf.WriteByte('\n')
		}
	}
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, buf.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	return p, seqs
}

// checkNoLocalPaths fails if the manifest at path holds an absolute path or
// any of the given local paths.
func checkNoLocalPaths(t *testing.T, path string, local ...string) {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, l := range append(local, os.TempDir()) {
		if l != "" && bytes.Contains(data, []byte(l)) {
			t.Errorf("manifest contains the local path %s", l)
		}
	}
	var doc any
	if err := json.Unmarshal(data, &doc); err != nil {
		t.Fatal(err)
	}
	var walk func(string, any)
	walk = func(at string, v any) {
		switch v := v.(type) {
		case map[string]any:
			for k, x := range v {
				walk(at+"."+k, x)
			}
		case []any:
			for _, x := range v {
				walk(at+"[]", x)
			}
		case string:
			if filepath.IsAbs(v) {
				t.Errorf("manifest %s is an absolute path: %s", at, v)
			}
		}
	}
	walk("", doc)
}

// checkInput fails unless in records the file at path.
func checkInput(t *testing.T, what string, in *manifest.Input, path string) {
	t.Helper()
	n, sum, err := manifest.HashFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if want := (manifest.Input{Name: filepath.Base(path), Bytes: n, SHA256: sum}); in == nil || *in != want {
		t.Errorf("manifest %s = %+v, want %+v", what, in, want)
	}
}

// A keyed run sampled from a reference region with planted variants is
// attributed only with its reference (and VCF): the manifest records the
// FASTA by hash, verify uses the file whose sha256 matches, and refuses
// without it or with another FASTA. --candidates based on the manifest read
// the same files.
func TestVerifyRegeneration_Reference(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := regenKeyFile(t, regenKey)
	inDir := t.TempDir()
	ref, seqs := writeFASTA(t, inDir, "ref.fa", 7)
	other, _ := writeFASTA(t, inDir, "other.fa", 70)
	alt := map[byte]string{'A': "C", 'C': "G", 'G': "T", 'T': "A"}
	vcf := filepath.Join(inDir, "calls.vcf")
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr2\t5000\t.\t" + string(seqs["chr2"][4999]) + "\t" + alt[seqs["chr2"][4999]] + "\t.\t.\tAF=0.5\n"
	if err := os.WriteFile(vcf, []byte(body), 0o600); err != nil {
		t.Fatal(err)
	}
	c := illuminaRun(t.TempDir())
	c.ReferencePath, c.ReferenceRegion, c.VariantsPath = ref, "chr2:1001-15000", vcf
	files, man := keyedRun(t, c, keyFile, "ACME", "study_r", 3)

	checkNoLocalPaths(t, man, inDir, ref, vcf, c.OutputDir, keyFile)
	m, err := manifest.Read(man)
	if err != nil {
		t.Fatal(err)
	}
	checkInput(t, "reference", m.Reference, ref)
	checkInput(t, "variants", m.Variants, vcf)
	if m.Config.ReferencePath != "" || m.Config.ReferenceRegion != "chr2:1001-15000" {
		t.Errorf("manifest config reference_path %q region %q", m.Config.ReferencePath, m.Config.ReferenceRegion)
	}

	ok := VerifyOptions{KeyFile: keyFile, Manifests: []string{man}, References: []string{other, ref}, Variants: []string{vcf}}
	out, err := verifyOut(t, files[:2], ok)
	if err != nil || strings.Count(out, "✅") != 2 || !strings.Contains(out, "200 of 200 reads (200 exact") {
		t.Errorf("with the reference and VCF: %v\n%s", err, out)
	}

	cands := filepath.Join(t.TempDir(), "cands.tsv")
	if err := os.WriteFile(cands, []byte("customer\tdataset\tsalt\nBETA\tstudy_r\t3\nACME\tstudy_r\t3\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	withCands := ok
	withCands.Candidates = cands
	if out, err := verifyOut(t, files[:1], withCands); err != nil || !strings.Contains(out, "✅") || !strings.Contains(out, "matches customer=ACME") {
		t.Errorf("--candidates on the manifest: %v\n%s", err, out)
	}
	withCands.References = nil
	if _, err := verifyOut(t, files[:1], withCands); err == nil || !strings.Contains(err.Error(), "--reference") {
		t.Errorf("--candidates without the reference: err = %v", err)
	}

	for name, tc := range map[string]struct {
		refs, vcfs []string
		want       []string
	}{
		"no reference":      {nil, []string{vcf}, []string{"ref.fa", "--reference", m.Reference.SHA256}},
		"another reference": {[]string{other}, []string{vcf}, []string{"ref.fa", "--reference"}},
		"no VCF":            {[]string{ref}, nil, []string{"calls.vcf", "--variants"}},
	} {
		o := ok
		o.References, o.Variants = tc.refs, tc.vcfs
		out, err := verifyOut(t, files[:1], o)
		if err == nil {
			t.Errorf("%s: verified\n%s", name, out)
			continue
		}
		for _, w := range tc.want {
			if !strings.Contains(err.Error(), w) {
				t.Errorf("%s: err = %v, want it to name %q", name, err, w)
			}
		}
	}
	if _, err := verifyOut(t, files[:1], VerifyOptions{KeyFile: keyFile, Manifests: []string{man}, References: []string{filepath.Join(inDir, "missing.fa")}}); err == nil || !strings.Contains(err.Error(), "--reference") {
		t.Errorf("missing --reference file: err = %v", err)
	}
}

// A manifest of phantom-science 0.14.0 and earlier records reference_path
// and no reference record: verify uses the file given with that name, or the
// recorded path, and warns that it cannot be checked.
func TestVerifyRegeneration_LegacyReferencePath(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := regenKeyFile(t, regenKey)
	inDir := t.TempDir()
	ref, _ := writeFASTA(t, inDir, "ref.fa", 7)
	c := illuminaRun(t.TempDir())
	c.ReferencePath = ref
	files, man := keyedRun(t, c, keyFile, "ACME", "study_l", 4)
	m, err := manifest.Read(man)
	if err != nil {
		t.Fatal(err)
	}
	legacy := func(path string) string {
		old := *m
		old.Reference = nil
		old.Config.ReferencePath = path
		data, err := old.Marshal()
		if err != nil {
			t.Fatal(err)
		}
		p := filepath.Join(t.TempDir(), manifest.FileName)
		if err := os.WriteFile(p, data, 0o600); err != nil {
			t.Fatal(err)
		}
		return p
	}
	copied := filepath.Join(t.TempDir(), "ref.fa")
	if err := os.WriteFile(copied, mustReadFile(t, ref), 0o600); err != nil {
		t.Fatal(err)
	}
	gone := filepath.Join(t.TempDir(), "moved", "ref.fa")
	for name, tc := range map[string]struct {
		recorded string
		refs     []string
		use      string
	}{
		"recorded path":         {ref, nil, ref},
		"same name given":       {gone, []string{ref, copied}, ref},
		"only file given":       {gone, []string{copied}, copied},
		"recorded path missing": {gone, nil, ""},
	} {
		out, err := verifyOut(t, files[:1], VerifyOptions{KeyFile: keyFile, Manifests: []string{legacy(tc.recorded)}, References: tc.refs})
		if tc.use == "" {
			if err == nil || !strings.Contains(err.Error(), "--reference") || !strings.Contains(err.Error(), "0.14.0 or earlier") {
				t.Errorf("%s: err = %v\n%s", name, err, out)
			}
			continue
		}
		if err != nil || !strings.Contains(out, "✅") || !strings.Contains(out, "[warn]") || !strings.Contains(out, "using "+tc.use+",") {
			t.Errorf("%s: %v\n%s", name, err, out)
		}
	}
}

func mustReadFile(t *testing.T, path string) []byte {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return data
}

// The HMM model and the empirical read-length histogram are recorded by
// hash too (no hmm_model_path, model_dir or empirical_path in the
// manifest) and verify needs them: --hmm-model, --read-length-histogram.
func TestVerifyRegeneration_ModelAndHistogram(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := regenKeyFile(t, regenKey)
	inDir := t.TempDir()
	model := authenticity.NewAdaptiveHMMModel(3)
	if err := model.TrainFromSequence(randomBases(5, 5000)); err != nil {
		t.Fatal(err)
	}
	if err := os.MkdirAll(filepath.Join(inDir, "models"), 0o750); err != nil {
		t.Fatal(err)
	}
	modelPath := filepath.Join(inDir, "models", "YRI.hmm.json")
	if err := model.Save(modelPath); err != nil {
		t.Fatal(err)
	}
	hist := filepath.Join(inDir, "lengths.tsv")
	if err := os.WriteFile(hist, mustReadFile(t, filepath.Join("..", "genomics", "readlen", "testdata", "cdna_hist.tsv")), 0o600); err != nil {
		t.Fatal(err)
	}
	other := filepath.Join(inDir, "other.tsv")
	if err := os.WriteFile(other, []byte("length\tcount\n300\t5\n900\t5\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	c := illuminaRun(t.TempDir())
	c.Population, c.ModelDir = "YRI", filepath.Join(inDir, "models")
	mFiles, mMan := keyedRun(t, c, keyFile, "ACME", "study_m", 8)
	o := ontRun(t.TempDir())
	o.ReadLengthDistribution = types.ReadLengthSpec{Dist: types.ReadLengthEmpirical, EmpiricalPath: hist}
	hFiles, hMan := keyedRun(t, o, keyFile, "ACME", "study_h", 9)

	for _, tc := range []struct {
		man   string
		files []string
		ok    VerifyOptions
		bad   []VerifyOptions
		flag  string
	}{
		{mMan, mFiles, VerifyOptions{HMMModels: []string{modelPath}}, []VerifyOptions{{}, {HMMModels: []string{hist}}}, "--hmm-model"},
		{hMan, hFiles, VerifyOptions{ReadLengths: []string{other, hist}}, []VerifyOptions{{}, {ReadLengths: []string{other}}}, "--read-length-histogram"},
	} {
		checkNoLocalPaths(t, tc.man, inDir, modelPath, hist)
		m, err := manifest.Read(tc.man)
		if err != nil {
			t.Fatal(err)
		}
		if m.Config.HMMModelPath != "" || m.Config.ModelDir != "" || m.Config.ReadLengthDistribution.EmpiricalPath != "" {
			t.Errorf("%s: manifest config keeps a path: %+v", tc.flag, m.Config)
		}
		if tc.flag == "--hmm-model" {
			checkInput(t, "hmm_model", m.HMMModel, modelPath)
		} else {
			checkInput(t, "read_length_histogram", m.ReadLengths, hist)
		}
		opts := tc.ok
		opts.KeyFile, opts.Manifests = keyFile, []string{tc.man}
		if out, err := verifyOut(t, tc.files[:1], opts); err != nil || !strings.Contains(out, "✅") {
			t.Errorf("%s: %v\n%s", tc.flag, err, out)
		}
		for _, bad := range tc.bad {
			bad.KeyFile, bad.Manifests = keyFile, []string{tc.man}
			if _, err := verifyOut(t, tc.files[:1], bad); err == nil || !strings.Contains(err.Error(), tc.flag) {
				t.Errorf("%s without the file: err = %v", tc.flag, err)
			}
		}
	}
}
