//go:build !release

package main

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"io/fs"
	"math/rand/v2"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"testing"
)

// The cross-architecture golden test pins the SHA-256 of the output of one
// small run of every generation mode. CI runs it on amd64 and on arm64, so
// any difference between CPU architectures (docs/CHARTER.md, principle 4)
// fails CI. The fixtures are built with integer arithmetic only, so they are
// the same on every machine. Set PHANTOM_SCIENCE_BIN to test a prebuilt
// binary instead of building one.
//
// When a change to the generators is intended to change output, re-pin by
// running with -v and copying the logged hashes; the change must then pass
// on both CI architectures with the same pins.

// goldenWant maps each case to the SHA-256 of its output. Pinned on
// darwin/arm64 at the change that moved the output path to internal/detmath
// (#33); the same values were produced there by v0.10.1.
var goldenWant = map[string]string{
	"illumina binned paired":               "20d5ea5a8d5dbf6f8ac06e550d0e96542672d3b049b7ea9d3e23f5ad32251019",
	"illumina reference paired":            "68ee5f9a8730e7bfab753839d5e0225af7dcab4860c5a8fc8193d53a3a7b6baf",
	"illumina lognormal inserts and reads": "cf08b7f96f5bd71ad86be9db2f731e9436514a9bcb9572c0fead5e435a136150",
	"ont hac lognormal":                    "be4106e6e7c7124f5bf7b9cb241f322e366c2c7b05f9c822828aa9713186e962",
	"ont fast gamma":                       "d953196a6882cd0a44f4dcf8e6400204d178fe4857a0cedb85b1a228fda3aa31",
	"signal dna":                           "3ba4105db252110bebf638a3bc18fc70c77b274c5ea1d3171a48f52c7cf1729b",
	"signal rna":                           "8b887aa607e7122281438c6b65e6103b90272cb3df6e343d2ef6bc65e75585c5",
	"transcriptome":                        "88585ae93756fa5a87a407761187024afc757fff36d0688f010d6e3af6bbe840",
	"scenario trio":                        "e1fc2b5f3b7c53d4585e9a91f4c4f38f4d7c49ed5711e9c84428da08f686be17",
	"hmm train and generate":               "324aee7b5f6caf72a035fd9a89b0a7bbf8398b91e343d0fc076f78ec5cfdd537",
}

func goldenBinary(t *testing.T) string {
	t.Helper()
	if bin := os.Getenv("PHANTOM_SCIENCE_BIN"); bin != "" {
		return bin
	}
	bin := filepath.Join(t.TempDir(), "phantom-science")
	if runtime.GOOS == "windows" {
		bin += ".exe"
	}
	build := exec.Command("go", "build", "-o", bin, ".")
	if out, err := build.CombinedOutput(); err != nil {
		t.Fatalf("go build: %v\n%s", err, out)
	}
	return bin
}

// goldenFixtures writes the inputs shared by the cases into dir.
type goldenFixtures struct {
	dir, seqConfig, pairedConfig, reference, poreModel, vcf, spec string
}

func writeGoldenFixtures(t *testing.T) goldenFixtures {
	t.Helper()
	dir := t.TempDir()
	write := func(name, content string) string {
		p := filepath.Join(dir, name)
		if err := os.WriteFile(p, []byte(content), 0o600); err != nil {
			t.Fatal(err)
		}
		return p
	}
	f := goldenFixtures{dir: dir}
	f.seqConfig = write("seq.json", `{"sample_name":"G","lane":1,"read_length":150,"reads_per_batch":500,"total_batches":1}`)
	f.pairedConfig = write("paired.json", `{"sample_name":"G","lane":1,"read_length":150,"reads_per_batch":500,"total_batches":1,"paired_end":true}`)

	// Reference: two 20 kb contigs of uniform random bases.
	r := rand.New(rand.NewPCG(33, 33))
	var ref strings.Builder
	for _, name := range []string{"chr1", "chr2"} {
		ref.WriteString(">" + name + "\n")
		for i := 0; i < 20000; i++ {
			ref.WriteByte("ACGT"[r.IntN(4)])
			if i%70 == 69 {
				ref.WriteByte('\n')
			}
		}
		ref.WriteByte('\n')
	}
	f.reference = write("ref.fa", ref.String())

	// Pore model: a synthetic k=5 table in pA, levels and spreads from
	// integer arithmetic, printed exactly.
	var pm strings.Builder
	pm.WriteString("#k\t5\nkmer\tlevel_mean\tlevel_stdv\n")
	for i := 0; i < 1<<10; i++ {
		b := make([]byte, 5)
		for j, idx := 4, i; j >= 0; j, idx = j-1, idx>>2 {
			b[j] = "ACGT"[idx&3]
		}
		lv := 450 + (i*617)%801 // tenths of a pA: 45.0–125.0
		fmt.Fprintf(&pm, "%s\t%d.%d\t1.%d\n", b, lv/10, lv%10, 3*(i%4))
	}
	f.poreModel = write("pore_5mer.tsv", pm.String())

	// VCF for HMM training.
	var vcf strings.Builder
	vcf.WriteString("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\n")
	refs := []string{"ATCG", "GCTA", "TATC", "CGAT", "AATG", "CCAT", "GGCA", "TTAC"}
	alts := []string{"ATCC", "GCTG", "TATA", "CGAC", "AATT", "CCAG", "GGCG", "TTAA"}
	for i := 0; i < 300; i++ {
		fmt.Fprintf(&vcf, "chr1\t%d\t.\t%s\t%s\n", (i+1)*10, refs[i%8], alts[(i/3)%8])
	}
	f.vcf = write("train.vcf", vcf.String())

	f.spec = write("trio.yaml", fmt.Sprintf(`phenotype:
  planted_variant: {chrom: chr1, pos: 2001, ref: C, alt: T, zygosity: heterozygous}
trio: {generate: true, de_novo: true}
sequencer: {platform: NovaSeq6000, coverage: 8, read_length: 100}
output:
  proband_fastq: %[1]s/proband.fastq
  mother_fastq: %[1]s/mother.fastq
  father_fastq: %[1]s/father.fastq
  ground_truth: %[1]s/truth.json
`, filepath.Join(dir, "trio")))
	return f
}

// hashTree returns the SHA-256 over the relative paths and contents of the
// files under dir, in sorted order, skipping files whose base name is in
// skip.
func hashTree(t *testing.T, dir string, skip ...string) string {
	t.Helper()
	var files []string
	err := filepath.WalkDir(dir, func(p string, d fs.DirEntry, err error) error {
		if err != nil || d.IsDir() {
			return err
		}
		for _, s := range skip {
			if d.Name() == s {
				return nil
			}
		}
		files = append(files, p)
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	sort.Strings(files)
	h := sha256.New()
	for _, p := range files {
		rel, _ := filepath.Rel(dir, p)
		data, err := os.ReadFile(p)
		if err != nil {
			t.Fatal(err)
		}
		fmt.Fprintf(h, "%s\x00%d\x00", filepath.ToSlash(rel), len(data))
		h.Write(data)
	}
	return hex.EncodeToString(h.Sum(nil))
}

func TestCrossArchGolden(t *testing.T) {
	if testing.Short() {
		t.Skip("builds the binary and runs every generation mode")
	}
	bin := goldenBinary(t)
	f := writeGoldenFixtures(t)
	run := func(t *testing.T, args ...string) []byte {
		t.Helper()
		cmd := exec.Command(bin, args...)
		cmd.Dir = f.dir
		var stdout, stderr bytes.Buffer
		cmd.Stdout, cmd.Stderr = &stdout, &stderr
		if err := cmd.Run(); err != nil {
			t.Fatalf("phantom-science %s: %v\n%s%s", strings.Join(args, " "), err, stdout.String(), stderr.String())
		}
		return stdout.Bytes()
	}
	sum := func(b []byte) string { h := sha256.Sum256(b); return hex.EncodeToString(h[:]) }
	ontReads := filepath.Join(f.dir, "ont.fastq")
	model := filepath.Join(f.dir, "model.hmm.json")
	// The run manifest records timestamps and paths; everything else is
	// hashed.
	manifests := []string{"phantom-manifest.json", "truth.json"}

	cases := []struct {
		name string
		run  func(t *testing.T) string
	}{
		{"illumina binned paired", func(t *testing.T) string {
			return sum(run(t, "generate", "--config", f.pairedConfig, "--stdout", "--seed", "7", "--total-reads", "1000"))
		}},
		{"illumina reference paired", func(t *testing.T) string {
			return sum(run(t, "generate", "--config", f.pairedConfig, "--stdout", "--seed", "7", "--total-reads", "1000",
				"--reference", f.reference, "--region", "chr2:1001-19000"))
		}},
		{"illumina lognormal inserts and reads", func(t *testing.T) string {
			return sum(run(t, "generate", "--config", f.pairedConfig, "--stdout", "--seed", "7", "--total-reads", "500",
				"--read-length-dist", "lognormal", "--read-length-mean", "140", "--read-length-sd", "15"))
		}},
		{"ont hac lognormal", func(t *testing.T) string {
			out := run(t, "generate", "--config", f.seqConfig, "--platform", "ont", "--basecaller-tier", "hac",
				"--stdout", "--seed", "7", "--total-reads", "60")
			if err := os.WriteFile(ontReads, out, 0o600); err != nil {
				t.Fatal(err)
			}
			return sum(out)
		}},
		{"ont fast gamma", func(t *testing.T) string {
			return sum(run(t, "generate", "--config", f.seqConfig, "--platform", "ont", "--basecaller-tier", "fast",
				"--stdout", "--seed", "7", "--total-reads", "60",
				"--read-length-dist", "gamma", "--read-length-mean", "3000", "--read-length-sd", "2000"))
		}},
		{"signal dna", func(t *testing.T) string {
			out := filepath.Join(f.dir, "sig-dna")
			run(t, "signal", "--input", ontReads, "--pore-model", f.poreModel, "--out", out, "--seed", "7")
			return hashTree(t, out, manifests...)
		}},
		{"signal rna", func(t *testing.T) string {
			out := filepath.Join(f.dir, "sig-rna")
			run(t, "signal", "--input", ontReads, "--pore-model", f.poreModel, "--out", out, "--seed", "7",
				"--profile", "rna004_130bps_4khz")
			return hashTree(t, out, manifests...)
		}},
		{"transcriptome", func(t *testing.T) string {
			out := filepath.Join(f.dir, "tx")
			fx := filepath.Join("internal", "genomics", "transcriptome", "testdata")
			wd, _ := os.Getwd()
			run(t, "transcriptome", "--genome", filepath.Join(wd, fx, "tiny_genome.fna"), "--gtf", filepath.Join(wd, fx, "tiny.gtf"),
				"--reads", "300", "--out", out, "--seed", "7", "--novel-rate", "0.3")
			return hashTree(t, out, manifests...)
		}},
		{"scenario trio", func(t *testing.T) string {
			cfg := filepath.Join(f.dir, "trio-seq.json")
			if err := os.WriteFile(cfg, []byte(`{"sample_name":"T","lane":1,"read_length":100,"seed":7}`), 0o600); err != nil {
				t.Fatal(err)
			}
			if err := os.MkdirAll(filepath.Join(f.dir, "trio"), 0o755); err != nil {
				t.Fatal(err)
			}
			run(t, "scenario", "--spec", f.spec, "--config", cfg, "--population", "YRI")
			return hashTree(t, filepath.Join(f.dir, "trio"), manifests...)
		}},
		{"hmm train and generate", func(t *testing.T) string {
			run(t, "train", "--vcf", f.vcf, "--out", model)
			data, err := os.ReadFile(model)
			if err != nil {
				t.Fatal(err)
			}
			// The model envelope records when it was serialised.
			if i := bytes.Index(data, []byte(`"serialized_at"`)); i >= 0 {
				if j := bytes.IndexByte(data[i:], ','); j > 0 {
					data = append(data[:i:i], data[i+j:]...)
				}
			}
			reads := run(t, "generate", "--config", f.seqConfig, "--hmm-model-path", model, "--stdout", "--seed", "7", "--total-reads", "300")
			return sum(append(data, reads...))
		}},
	}
	for _, tc := range cases {
		got := tc.run(t)
		want, ok := goldenWant[tc.name]
		switch {
		case !ok:
			t.Errorf("%s: no pinned hash; got %s", tc.name, got)
		case got != want:
			t.Errorf("%s: output hash %s, want %s (on %s/%s)", tc.name, got, want, runtime.GOOS, runtime.GOARCH)
		default:
			t.Logf("%s: %s", tc.name, got)
		}
	}
}
