//go:build !release

package cmd

// Keyed-seed fingerprints of the scenario, signal and transcriptome
// commands (#31), end to end: a fingerprinted run is the unfingerprinted
// run with the derived seed (so its data and statistics are those of any
// run), its manifest records the fingerprint and no seed, key or path, and
// verify attributes its output to the right customer only.

import (
	"bytes"
	"fmt"
	"math/rand/v2"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/signalrun"
	"github.com/scttfrdmn/phantom-science/internal/genomics/slow5"
	"github.com/spf13/pflag"
)

const kindsKey = "kinds-verify-key-0123456789abcdef" // gitleaks:allow -- test fixture

// runCLI runs "phantom-science args..." and resets every flag afterwards.
func runCLI(t *testing.T, args ...string) (stdout, stderr string, err error) {
	t.Helper()
	var out, errb bytes.Buffer
	rootCmd.SetOut(&out)
	rootCmd.SetErr(&errb)
	rootCmd.SetArgs(args)
	defer func() {
		rootCmd.SetArgs([]string{})
		rootCmd.SetOut(nil)
		rootCmd.SetErr(nil)
		reset := func(f *pflag.Flag) {
			if sv, ok := f.Value.(pflag.SliceValue); ok {
				_ = sv.Replace(nil)
			} else {
				_ = f.Value.Set(f.DefValue)
			}
			f.Changed = false
		}
		rootCmd.PersistentFlags().VisitAll(reset)
		for _, c := range rootCmd.Commands() {
			c.Flags().VisitAll(reset)
		}
	}()
	err = rootCmd.Execute()
	return out.String(), errb.String(), err
}

func writeKey(t *testing.T, key string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "fp.key")
	if err := os.WriteFile(p, []byte(key+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

// fpConfig writes a sequencer config enabling the keyed-seed fingerprint
// for customer and dataset; extra JSON fields are added verbatim.
func fpConfig(t *testing.T, keyFile, customer, dataset, extra string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "fp.json")
	body := fmt.Sprintf(`{%s"fingerprinting":{"enable":true,"key_file":%q},"enterprise":{"customer_id":%q},"dataset_id":%q}`,
		extra, keyFile, customer, dataset)
	if err := os.WriteFile(p, []byte(body), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

func derived(t *testing.T, customer, dataset string, salt uint64) uint64 {
	t.Helper()
	s, err := fingerprinting.DeriveSeed([]byte(kindsKey), customer, dataset, salt)
	if err != nil {
		t.Fatal(err)
	}
	return s
}

// checkPrivateManifest checks that the manifest at path records the keyed
// fingerprint (customer, dataset, salt) and neither a seed nor any of the
// forbidden strings (key, derived seed, local paths).
func checkPrivateManifest(t *testing.T, path, customer, dataset string, salt uint64, forbidden ...string) *manifest.Manifest {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	m, err := manifest.Read(path)
	if err != nil {
		t.Fatal(err)
	}
	want := manifest.Fingerprint{Scheme: fingerprinting.KeyedSeedScheme, Customer: customer, Dataset: dataset, Salt: salt}
	if m.Fingerprint == nil || *m.Fingerprint != want {
		t.Errorf("%s: fingerprint %+v, want %+v", path, m.Fingerprint, want)
	}
	if m.Seed != 0 {
		t.Errorf("%s records seed %d", path, m.Seed)
	}
	for _, fc := range m.FlowCells {
		if fc.Seed != 0 {
			t.Errorf("%s: flow cell %d records seed %d", path, fc.Index, fc.Seed)
		}
	}
	if bytes.Contains(data, []byte(`"seed": `)) && m.Config.Seed != salt {
		t.Errorf("%s: config seed %d, want the salt %d", path, m.Config.Seed, salt)
	}
	for _, f := range append(forbidden, kindsKey) {
		if f != "" && bytes.Contains(data, []byte(f)) {
			t.Errorf("%s contains %q", path, f)
		}
	}
	return m
}

func verifyKinds(t *testing.T, paths []string, opts cli.VerifyOptions) (string, error) {
	t.Helper()
	var out bytes.Buffer
	err := cli.VerifyFiles(&out, paths, opts)
	return out.String(), err
}

// sameFiles fails unless dirs a and b hold the same files with the same
// bytes, apart from those named in skip.
func sameFiles(t *testing.T, a, b string, skip ...string) {
	t.Helper()
	fa, fb := readDirFiles(t, a), readDirFiles(t, b)
	for _, s := range skip {
		delete(fa, s)
		delete(fb, s)
	}
	if len(fa) != len(fb) || len(fa) == 0 {
		t.Fatalf("%s has %d files, %s %d", a, len(fa), b, len(fb))
	}
	for name, data := range fa {
		if !bytes.Equal(data, fb[name]) {
			t.Errorf("%s differs between %s and %s", name, a, b)
		}
	}
}

func TestSignalFingerprint(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := writeKey(t, kindsKey)
	inputs := generateONTReads(t)
	model := writeTestPoreModel(t)
	cfgA := fpConfig(t, keyFile, "ACME", "sig_1", "")
	cfgB := fpConfig(t, keyFile, "BETA", "sig_1", "")
	base := t.TempDir()
	run := func(name string, args ...string) string {
		t.Helper()
		out := filepath.Join(base, name)
		a := []string{"signal", "--input", inputs[0], "--input", inputs[1], "--pore-model", model, "--out", out, "--compression", "none"}
		if _, stderr, err := runCLI(t, append(a, args...)...); err != nil {
			t.Fatalf("signal %v: %v\n%s", args, err, stderr)
		}
		return out
	}
	a := run("a", "--config", cfgA, "--seed", "9", "--workers", "1")
	a4 := run("a4", "--config", cfgA, "--seed", "9", "--workers", "4")
	plain := run("plain", "--seed", strconv.FormatUint(derived(t, "ACME", "sig_1", 9), 10))
	b := run("b", "--config", cfgB, "--seed", "9")

	// The fingerprinted run is the unfingerprinted run with the derived
	// seed, for any number of workers; only the manifests differ.
	sameFiles(t, a, a4)
	sameFiles(t, a, plain, manifest.FileName)
	manA := filepath.Join(a, manifest.FileName)
	checkPrivateManifest(t, manA, "ACME", "sig_1", 9,
		strconv.FormatUint(derived(t, "ACME", "sig_1", 9), 10), keyFile, filepath.Dir(keyFile), filepath.Dir(inputs[0]), filepath.Dir(model), base)
	pm, err := manifest.Read(filepath.Join(plain, manifest.FileName))
	if err != nil || pm.Fingerprint != nil || pm.Seed != derived(t, "ACME", "sig_1", 9) {
		t.Errorf("unfingerprinted manifest: %v %+v", err, pm)
	}

	blowA, _ := filepathGlob(t, a, "*.blow5")
	blowB, _ := filepathGlob(t, b, "*.blow5")
	_, _, recs := readBLOW5Dir(t, a)
	both := cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{filepath.Join(b, manifest.FileName), manA},
		SignalInputs: inputs, PoreModels: []string{model}}

	// Renamed reads, trimmed signal, repacked as SLOW5 ASCII.
	transformed := filepath.Join(t.TempDir(), "renamed.slow5")
	{
		f, err := os.Open(blowA)
		if err != nil {
			t.Fatal(err)
		}
		h, rs, err := slow5.ReadAll(f)
		f.Close()
		if err != nil {
			t.Fatal(err)
		}
		var buf bytes.Buffer
		w, err := slow5.NewSLOW5Writer(&buf, h)
		if err != nil {
			t.Fatal(err)
		}
		for i, r := range rs {
			r.ReadID = signalrun.FormatUUID(signalrun.UUIDv5(signalrun.Namespace, fmt.Sprintf("renamed-%d", i)))
			r.RawSignal = r.RawSignal[200 : len(r.RawSignal)-100]
			if err := w.Write(r); err != nil {
				t.Fatal(err)
			}
		}
		if err := w.Close(); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(transformed, buf.Bytes(), 0o600); err != nil {
			t.Fatal(err)
		}
	}
	n := len(recs)
	for name, tc := range map[string]struct {
		path string
		want string
	}{
		"as written":                  {blowA, fmt.Sprintf("%d of %d reads by raw signal (%d exact, 0 by", n, n, n)},
		"renamed, trimmed, SLOW5":     {transformed, fmt.Sprintf("%d of %d reads by raw signal (0 exact, %d by", n, n, n)},
		"customer B's own BLOW5 file": {blowB, "customer=BETA"},
	} {
		t.Run(name, func(t *testing.T) {
			out, err := verifyKinds(t, []string{tc.path}, both)
			if err != nil || !strings.Contains(out, "✅") || !strings.Contains(out, tc.want) || strings.Contains(out, "also matches") {
				t.Errorf("%v\n%s", err, out)
			}
			if tc.path != blowB && !strings.Contains(out, "matches a/phantom-manifest.json customer=ACME") {
				t.Errorf("not attributed to ACME:\n%s", out)
			}
			if !strings.Contains(out, "decoy run 0 (0 exact, 0 whole)") {
				t.Errorf("want a decoy run with no matches:\n%s", out)
			}
		})
	}

	cands := filepath.Join(t.TempDir(), "cands.tsv")
	if err := os.WriteFile(cands, []byte("BETA\tsig_1\t9\nACME\tsig_1\t10\nACME\tsig_2\t9\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	for name, opts := range map[string]cli.VerifyOptions{
		"customer B only":       {KeyFile: keyFile, Manifests: []string{filepath.Join(b, manifest.FileName)}, SignalInputs: inputs, PoreModels: []string{model}},
		"wrong key":             {KeyFile: writeKey(t, "another-key-0123456789abcdef"), Manifests: []string{manA}, SignalInputs: inputs, PoreModels: []string{model}},
		"wrong salt or dataset": {KeyFile: keyFile, Manifests: []string{manA}, Candidates: cands, SignalInputs: inputs, PoreModels: []string{model}},
	} {
		t.Run(name, func(t *testing.T) {
			out, err := verifyKinds(t, []string{transformed}, opts)
			if err == nil || strings.Contains(out, "✅") || !strings.Contains(out, " 0 of ") {
				t.Errorf("matched: %v\n%s", err, out)
			}
		})
	}
	for name, tc := range map[string]struct {
		opts cli.VerifyOptions
		want string
	}{
		"no input reads": {cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{manA}, PoreModels: []string{model}}, "pass that read file with --signal-input"},
		"no pore model":  {cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{manA}, SignalInputs: inputs}, "pass that pore model with --pore-model"},
		"other reads":    {cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{manA}, SignalInputs: inputs[:1], PoreModels: []string{model}}, "--signal-input"},
	} {
		t.Run(name, func(t *testing.T) {
			if out, err := verifyKinds(t, []string{blowA}, tc.opts); err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Errorf("err = %v, want %q\n%s", err, tc.want, out)
			}
		})
	}
	// A FASTQ suspect needs a candidate that wrote reads.
	if _, err := verifyKinds(t, inputs[:1], both); err == nil || !strings.Contains(err.Error(), "every candidate is a signal run") {
		t.Errorf("FASTQ against signal candidates: %v", err)
	}

	// Without --seed the salt is drawn, logged and recorded; it reproduces
	// the run.
	out := filepath.Join(base, "drawn")
	_, stderr, err := runCLI(t, "signal", "--input", inputs[0], "--pore-model", model, "--out", out, "--config", cfgA)
	if err != nil {
		t.Fatal(err)
	}
	mt := regexp.MustCompile(`fingerprint salt (\d+)`).FindStringSubmatch(stderr)
	if mt == nil {
		t.Fatalf("no salt logged:\n%s", stderr)
	}
	salt, _ := strconv.ParseUint(mt[1], 10, 64)
	checkPrivateManifest(t, filepath.Join(out, manifest.FileName), "ACME", "sig_1", salt, strconv.FormatUint(derived(t, "ACME", "sig_1", salt), 10))
	drawnFiles, _ := filepathGlob(t, out, "*.blow5")
	if res, err := verifyKinds(t, []string{drawnFiles}, cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{filepath.Join(out, manifest.FileName)}, SignalInputs: inputs, PoreModels: []string{model}}); err != nil {
		t.Errorf("drawn salt: %v\n%s", err, res)
	}
}

// filepathGlob returns the one file of dir matching pattern.
func filepathGlob(t *testing.T, dir, pattern string) (string, []string) {
	t.Helper()
	m, err := filepath.Glob(filepath.Join(dir, pattern))
	if err != nil || len(m) == 0 {
		t.Fatalf("%s/%s: %v %v", dir, pattern, m, err)
	}
	return m[0], m
}

func TestTranscriptomeFingerprint(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := writeKey(t, kindsKey)
	genome, gtf := filepath.Join(txFixtures, "tiny_genome.fna"), filepath.Join(txFixtures, "tiny.gtf")
	base := t.TempDir()
	run := func(name string, args ...string) string {
		t.Helper()
		out := filepath.Join(base, name)
		a := []string{"transcriptome", "--genome", genome, "--gtf", gtf, "--reads", "300", "--out", out, "--novel-rate", "0.3"}
		if _, stderr, err := runCLI(t, append(a, args...)...); err != nil {
			t.Fatalf("transcriptome %v: %v\n%s", args, err, stderr)
		}
		return out
	}
	a := run("a", "--config", fpConfig(t, keyFile, "ACME", "", ""), "--dataset-id", "tx_1", "--seed", "4")
	plain := run("plain", "--seed", strconv.FormatUint(derived(t, "ACME", "tx_1", 4), 10))
	b := run("b", "--config", fpConfig(t, keyFile, "BETA", "tx_1", ""), "--seed", "4")
	sameFiles(t, a, plain, manifest.FileName)
	manA := filepath.Join(a, manifest.FileName)
	m := checkPrivateManifest(t, manA, "ACME", "tx_1", 4,
		strconv.FormatUint(derived(t, "ACME", "tx_1", 4), 10), keyFile, filepath.Dir(keyFile), txFixtures, base)
	if m.Transcriptome == nil || m.Transcriptome.Genome == nil || m.Transcriptome.GTF == nil || m.Transcriptome.Transcripts != nil || m.Totals.Reads != 300 {
		t.Fatalf("manifest %+v", m)
	}
	if err := m.VerifyDir(a); err != nil {
		t.Errorf("outputs: %v", err)
	}

	reads := filepath.Join(a, "reads.fastq")
	recs := readFQ(t, reads)
	var renamed, revcomped bytes.Buffer
	for i, r := range recs {
		fmt.Fprintf(&renamed, "@read%d\n%s\n+\n%s\n", i, r[1], strings.Repeat("I", len(r[1])))
		fmt.Fprintf(&revcomped, "@read%d\n%s\n+\n%s\n", i, revcomp(r[1]), r[2])
	}
	write := func(name string, data []byte) string {
		p := filepath.Join(t.TempDir(), name)
		if err := os.WriteFile(p, data, 0o600); err != nil {
			t.Fatal(err)
		}
		return p
	}
	both := cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{filepath.Join(b, manifest.FileName), manA}, Genomes: []string{genome}, GTFs: []string{gtf}}
	for name, tc := range map[string]struct{ path, want string }{
		"as written":                   {reads, "300 of 300 reads (300 exact"},
		"renamed, qualities rewritten": {write("renamed.fastq", renamed.Bytes()), "300 of 300 reads (300 exact"},
		"reverse complemented":         {write("rc.fastq", revcomped.Bytes()), "300 of 300 reads (300 exact"},
	} {
		t.Run(name, func(t *testing.T) {
			out, err := verifyKinds(t, []string{tc.path}, both)
			if err != nil || !strings.Contains(out, "✅") || !strings.Contains(out, "matches a/phantom-manifest.json customer=ACME") ||
				!strings.Contains(out, tc.want) || strings.Contains(out, "also matches") || !strings.Contains(out, "decoy run") {
				t.Errorf("%v\n%s", err, out)
			}
		})
	}
	// Reads sampled from the same transcripts share their k-mers with any
	// run, so the wrong customer's run contains them; the decoy run shows
	// that, and they are not attributed to it.
	out, err := verifyKinds(t, []string{reads}, cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{filepath.Join(b, manifest.FileName)}, Genomes: []string{genome}, GTFs: []string{gtf}})
	if err == nil || strings.Contains(out, "✅") {
		t.Errorf("customer B: %v\n%s", err, out)
	}
	if _, err := verifyKinds(t, []string{reads}, cli.VerifyOptions{KeyFile: writeKey(t, "another-key-0123456789abcdef"), Manifests: []string{manA}, Genomes: []string{genome}, GTFs: []string{gtf}}); err == nil {
		t.Error("wrong key matched")
	}
	if _, err := verifyKinds(t, []string{reads}, cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{manA}, Genomes: []string{genome}}); err == nil || !strings.Contains(err.Error(), "pass that GTF with --gtf") {
		t.Errorf("missing GTF: %v", err)
	}
}

// readFQ returns the header, sequence and quality of every record of a
// plain FASTQ file.
func readFQ(t *testing.T, path string) [][3]string {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSuffix(string(data), "\n"), "\n")
	var out [][3]string
	for i := 0; i+3 < len(lines); i += 4 {
		out = append(out, [3]string{lines[i], lines[i+1], lines[i+3]})
	}
	return out
}

// writeScenarioSpec writes a trio scenario spec (or, with ref set, a
// proband-only one whose planted REF matches chr1 of the reference) whose
// outputs go to dir.
func writeScenarioSpec(t *testing.T, dir string, refBase byte) string {
	t.Helper()
	if err := os.MkdirAll(dir, 0o750); err != nil {
		t.Fatal(err)
	}
	var spec string
	if refBase == 0 {
		spec = fmt.Sprintf(`phenotype:
  planted_variant: {chrom: chr1, pos: 2001, ref: C, alt: T, zygosity: heterozygous}
trio: {generate: true, de_novo: true}
sequencer: {platform: NovaSeq6000, coverage: 6, read_length: 100}
output:
  proband_fastq: %[1]s/proband.fastq
  mother_fastq: %[1]s/mother.fastq
  father_fastq: %[1]s/father.fastq
  ground_truth: %[1]s/truth.json
`, dir)
	} else {
		alt := "ACGT"[(strings.IndexByte("ACGT", refBase)+1)%4]
		spec = fmt.Sprintf(`phenotype:
  planted_variant: {chrom: chr1, pos: 4001, ref: %c, alt: %c, zygosity: heterozygous}
sequencer: {platform: NovaSeq6000, coverage: 10, read_length: 100, paired_end: true}
output: {proband_fastq: %[3]s/p_R1.fastq, proband_fastq_r2: %[3]s/p_R2.fastq, ground_truth: %[3]s/truth.json}
`, refBase, alt, dir)
	}
	p := filepath.Join(t.TempDir(), "spec.yaml")
	if err := os.WriteFile(p, []byte(spec), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestScenarioFingerprint(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := writeKey(t, kindsKey)
	base := t.TempDir()
	const extra = `"sample_name":"T","lane":1,"read_length":100,"seed":7,`
	run := func(name, cfg string, args ...string) string {
		t.Helper()
		dir := filepath.Join(base, name)
		a := []string{"scenario", "--spec", writeScenarioSpec(t, dir, 0), "--config", cfg, "--population", "YRI"}
		if _, stderr, err := runCLI(t, append(a, args...)...); err != nil {
			t.Fatalf("scenario %s: %v\n%s", name, err, stderr)
		}
		return dir
	}
	plainCfg := filepath.Join(t.TempDir(), "plain.json")
	if err := os.WriteFile(plainCfg, []byte(fmt.Sprintf(`{"sample_name":"T","lane":1,"read_length":100,"seed":%d}`, derived(t, "ACME", "trio_1", 7))), 0o600); err != nil {
		t.Fatal(err)
	}
	a := run("a", fpConfig(t, keyFile, "ACME", "trio_1", extra))
	plain := run("plain", plainCfg)
	b := run("b", fpConfig(t, keyFile, "BETA", "trio_1", extra))

	// Every member's reads and the ground truth are those of the
	// unfingerprinted scenario with the derived seed.
	sameFiles(t, a, plain, manifest.FileName)
	manA := filepath.Join(a, manifest.FileName)
	m := checkPrivateManifest(t, manA, "ACME", "trio_1", 7,
		strconv.FormatUint(derived(t, "ACME", "trio_1", 7), 10), keyFile, filepath.Dir(keyFile), base)
	if m.Scenario == nil || strings.Join(m.Scenario.Members, ",") != "proband,mother,father" || m.Totals.Files != 3 || len(m.Outputs) != 4 {
		t.Fatalf("manifest %+v", m)
	}
	if err := m.VerifyDir(a); err != nil {
		t.Errorf("outputs: %v", err)
	}
	pm, err := manifest.Read(filepath.Join(plain, manifest.FileName))
	if err != nil || pm.Fingerprint != nil || pm.Seed != derived(t, "ACME", "trio_1", 7) {
		t.Errorf("unfingerprinted manifest: %v %+v", err, pm)
	}

	both := cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{filepath.Join(b, manifest.FileName), manA}}
	for _, member := range []string{"proband", "mother", "father"} {
		recs := readFQ(t, filepath.Join(a, member+".fastq"))
		var renamed bytes.Buffer
		for i, r := range recs {
			fmt.Fprintf(&renamed, "@x%d\n%s\n+\n%s\n", i, r[1][5:len(r[1])-5], strings.Repeat("F", len(r[1])-10))
		}
		p := filepath.Join(t.TempDir(), member+".fastq")
		if err := os.WriteFile(p, renamed.Bytes(), 0o600); err != nil {
			t.Fatal(err)
		}
		for _, path := range []string{filepath.Join(a, member+".fastq"), p} {
			out, err := verifyKinds(t, []string{path}, both)
			if err != nil || !strings.Contains(out, "✅") || !strings.Contains(out, "matches a/phantom-manifest.json customer=ACME") || strings.Contains(out, "also matches") {
				t.Errorf("%s: %v\n%s", path, err, out)
			}
		}
	}
	proband := filepath.Join(a, "proband.fastq")
	if out, err := verifyKinds(t, []string{proband}, cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{filepath.Join(b, manifest.FileName)}}); err == nil || !strings.Contains(out, " 0 of ") {
		t.Errorf("customer B: %v\n%s", err, out)
	}
	if out, err := verifyKinds(t, []string{proband}, cli.VerifyOptions{KeyFile: writeKey(t, "another-key-0123456789abcdef"), Manifests: []string{manA}}); err == nil {
		t.Errorf("wrong key matched:\n%s", out)
	}

	// A reference scenario is regenerated only with that reference.
	ref, seqs := writeCmdReference(t)
	dir := filepath.Join(base, "ref")
	if _, stderr, err := runCLI(t, "scenario", "--spec", writeScenarioSpec(t, dir, seqs["chr1"][4000]), "--reference", ref,
		"--config", fpConfig(t, keyFile, "ACME", "ref_1", `"sample_name":"P","lane":1,"read_length":100,"seed":3,`)); err != nil {
		t.Fatalf("reference scenario: %v\n%s", err, stderr)
	}
	refMan := filepath.Join(dir, manifest.FileName)
	rm := checkPrivateManifest(t, refMan, "ACME", "ref_1", 3, ref, filepath.Dir(ref))
	if rm.Reference == nil || rm.Reference.Name != "ref.fa" {
		t.Fatalf("reference record %+v", rm.Reference)
	}
	r1 := filepath.Join(dir, "p_R1.fastq")
	if _, err := verifyKinds(t, []string{r1}, cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{refMan}}); err == nil || !strings.Contains(err.Error(), "pass that FASTA with --reference") {
		t.Errorf("without --reference: %v", err)
	}
	other, _ := writeOtherReference(t)
	if _, err := verifyKinds(t, []string{r1}, cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{refMan}, References: []string{other}}); err == nil || !strings.Contains(err.Error(), "--reference") {
		t.Errorf("with another reference: %v", err)
	}
	out, err := verifyKinds(t, []string{r1, filepath.Join(dir, "p_R2.fastq")}, cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{refMan}, References: []string{other, ref}})
	if err != nil || strings.Count(out, "✅") != 2 || !strings.Contains(out, "decoy run") {
		t.Errorf("with the reference: %v\n%s", err, out)
	}
}

// writeOtherReference writes a FASTA unlike writeCmdReference's.
func writeOtherReference(t *testing.T) (string, string) {
	t.Helper()
	r := rand.New(rand.NewPCG(99, 99))
	b := make([]byte, 8000)
	for i := range b {
		b[i] = "ACGT"[r.IntN(4)]
	}
	p := filepath.Join(t.TempDir(), "other.fa")
	if err := os.WriteFile(p, []byte(">chr1\n"+string(b)+"\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	return p, string(b)
}

// The commands without a sequencer configuration refuse fingerprint
// settings they cannot honour, before writing anything.
func TestFingerprintSettingsRefused(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := writeKey(t, kindsKey)
	genome, gtf := filepath.Join(txFixtures, "tiny_genome.fna"), filepath.Join(txFixtures, "tiny.gtf")
	for name, tc := range map[string]struct {
		cfg, want string
	}{
		"no dataset":  {fpConfig(t, keyFile, "ACME", "", ""), "needs a dataset ID"},
		"no customer": {fpConfig(t, keyFile, "", "d", ""), "needs a customer ID"},
		"no key":      {fpConfig(t, "", "ACME", "d", ""), "needs a key"},
	} {
		t.Run(name, func(t *testing.T) {
			out := filepath.Join(t.TempDir(), "tx")
			_, _, err := runCLI(t, "transcriptome", "--genome", genome, "--gtf", gtf, "--reads", "10", "--out", out, "--config", tc.cfg)
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Errorf("err = %v, want %q", err, tc.want)
			}
			if _, err := os.Stat(out); err == nil {
				t.Errorf("%s was written", out)
			}
		})
	}
	labels := filepath.Join(t.TempDir(), "labels.json")
	if err := os.WriteFile(labels, []byte(fmt.Sprintf(`{"fingerprinting":{"enable":true,"key_file":%q,"customer_fingerprint":true},"enterprise":{"customer_id":"A"},"dataset_id":"d"}`, keyFile)), 0o600); err != nil {
		t.Fatal(err)
	}
	inputs := generateONTReads(t)
	if _, _, err := runCLI(t, "signal", "--input", inputs[0], "--pore-model", writeTestPoreModel(t), "--out", t.TempDir(), "--config", labels); err == nil ||
		!strings.Contains(err.Error(), "supports only the \"keyed_seed\" fingerprint layer") {
		t.Errorf("labels: %v", err)
	}
}

// Trimmed cDNA reads are not exact, but every sampled k-mer of each, its
// sequencing errors included, is in the run that made it (whole), and in
// few reads of another run of the same transcripts: they are attributed by
// the whole-match test against the decoy run.
func TestTranscriptomeFingerprint_Trimmed(t *testing.T) {
	t.Setenv(fingerprinting.KeyEnvVar, "")
	keyFile := writeKey(t, kindsKey)
	r := rand.New(rand.NewPCG(8, 8))
	var tx bytes.Buffer
	for i := range 200 {
		b := make([]byte, 800+r.IntN(1200))
		for j := range b {
			b[j] = "ACGT"[r.IntN(4)]
		}
		fmt.Fprintf(&tx, ">TX%03d|GENE%03d\n%s\n", i, i/2, b)
	}
	txPath := filepath.Join(t.TempDir(), "tx.fa")
	if err := os.WriteFile(txPath, tx.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	base := t.TempDir()
	run := func(name, customer string) string {
		t.Helper()
		out := filepath.Join(base, name)
		if _, stderr, err := runCLI(t, "transcriptome", "--transcripts", txPath, "--reads", "2000", "--out", out, "--seed", "6",
			"--config", fpConfig(t, keyFile, customer, "tx_2", "")); err != nil {
			t.Fatalf("%v\n%s", err, stderr)
		}
		return out
	}
	a, b := run("a", "ACME"), run("b", "BETA")
	var trimmed bytes.Buffer
	for i, rec := range readFQ(t, filepath.Join(a, "reads.fastq")) {
		fmt.Fprintf(&trimmed, "@t%d\n%s\n+\n%s\n", i, rec[1][10:len(rec[1])-10], strings.Repeat("F", len(rec[1])-20))
	}
	p := filepath.Join(t.TempDir(), "trimmed.fastq")
	if err := os.WriteFile(p, trimmed.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	opts := cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{filepath.Join(b, manifest.FileName), filepath.Join(a, manifest.FileName)}, Transcripts: []string{txPath}}
	out, err := verifyKinds(t, []string{p}, opts)
	if err != nil || !strings.Contains(out, "✅") || !strings.Contains(out, "matches a/phantom-manifest.json customer=ACME") ||
		!strings.Contains(out, "2000 of 2000 reads (0 exact, 2000 by 31-mer containment); 2000 whole") || strings.Contains(out, "also matches") {
		t.Errorf("%v\n%s", err, out)
	}
	opts.Manifests = opts.Manifests[:1]
	if out, err := verifyKinds(t, []string{p}, opts); err == nil {
		t.Errorf("customer B matched:\n%s", out)
	}
}
