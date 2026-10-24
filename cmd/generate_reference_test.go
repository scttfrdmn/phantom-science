//go:build !release

package cmd

import (
	"bytes"
	"fmt"
	"io"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
)

// writeCmdReference writes a two-contig FASTA and returns its path and
// sequences.
func writeCmdReference(t *testing.T) (string, map[string]string) {
	t.Helper()
	r := rand.New(rand.NewPCG(9, 9))
	seqs := map[string]string{}
	var buf bytes.Buffer
	for _, name := range []string{"chr1", "chr2"} {
		b := make([]byte, 8000)
		for i := range b {
			b[i] = "ACGT"[r.IntN(4)]
		}
		seqs[name] = string(b)
		buf.WriteString(">" + name + "\n")
		for i := 0; i < len(b); i += 70 {
			buf.Write(b[i:min(i+70, len(b))])
			buf.WriteByte('\n')
		}
	}
	path := filepath.Join(t.TempDir(), "ref.fa")
	if err := os.WriteFile(path, buf.Bytes(), 0o600); err != nil {
		t.Fatal(err)
	}
	return path, seqs
}

// writePairedConfig writes a paired-end config, with extra JSON fields.
func writePairedConfig(t *testing.T, extra string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "seq.json")
	cfg := fmt.Sprintf(`{"sample_name":"REF","lane":1,"output_dir":%q,"read_length":100,"reads_per_batch":200,"total_batches":2,"paired_end":true%s}`,
		filepath.Join(dir, "out"), extra)
	if err := os.WriteFile(path, []byte(cfg), 0o600); err != nil {
		t.Fatal(err)
	}
	return path
}

func revcomp(s string) string {
	c := map[byte]byte{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	b := make([]byte, len(s))
	for i := range s {
		b[len(s)-1-i] = c[s[i]]
	}
	return string(b)
}

// generate --reference --region in file mode samples pairs from the region
// and writes truth tables listed in the manifest.
func TestGenerate_ReferenceRegionFiles(t *testing.T) {
	ref, seqs := writeCmdReference(t)
	out := filepath.Join(t.TempDir(), "out")
	setGenerateFlags(t, map[string]string{
		"config": writePairedConfig(t, ""), "reference": ref, "region": "chr2:2001-6000",
		"illumina-errors": "none", "seed": "5", "output": out,
	})
	if err := runGenerate(generateCmd, io.Discard); err != nil {
		t.Fatal(err)
	}
	m, err := manifest.Read(filepath.Join(out, manifest.FileName))
	if err != nil {
		t.Fatal(err)
	}
	if err := m.VerifyDir(out); err != nil {
		t.Fatal(err)
	}
	if m.Config.ReferencePath != ref || m.Config.ReferenceRegion != "chr2:2001-6000" {
		t.Errorf("manifest config reference %q region %q", m.Config.ReferencePath, m.Config.ReferenceRegion)
	}
	truths := 0
	for _, o := range m.Outputs {
		if o.Role != manifest.RoleTruth {
			continue
		}
		truths++
		data, err := os.ReadFile(filepath.Join(out, o.Path))
		if err != nil {
			t.Fatal(err)
		}
		r1, err := os.ReadFile(filepath.Join(out, strings.Replace(o.Path, ".truth.tsv", ".fastq", 1)))
		if err != nil {
			t.Fatal(err)
		}
		lines := strings.Split(strings.TrimSpace(string(data)), "\n")[1:]
		fq := strings.Split(strings.TrimSpace(string(r1)), "\n")
		if len(lines) != 200 || len(fq) != 4*200 {
			t.Fatalf("%s: %d rows, %d FASTQ lines", o.Path, len(lines), len(fq))
		}
		for i, line := range lines {
			f := strings.Split(line, "\t")
			start, _ := strconv.ParseInt(f[2], 10, 64)
			n, _ := strconv.Atoi(f[4])
			if f[1] != "chr2" || start < 2000 || start+int64(n) > 6000 {
				t.Fatalf("fragment %s outside chr2:2001-6000", line)
			}
			frag := seqs["chr2"][start : start+int64(n)]
			want := frag[:100]
			if f[3] == "-" {
				want = revcomp(frag[n-100:])
			}
			if got := fq[4*i+1]; got != want {
				t.Fatalf("%s row %d: R1 is not the reference fragment", o.Path, i)
			}
		}
	}
	if truths != 2 {
		t.Errorf("%d truth tables, want 2", truths)
	}
}

// With --reference, a VCF REF that does not match fails with --strict-ref;
// --reference with an HMM model is a configuration error.
func TestGenerate_ReferenceErrors(t *testing.T) {
	ref, seqs := writeCmdReference(t)
	vcf := filepath.Join(t.TempDir(), "v.vcf")
	next := map[string]string{"A": "C", "C": "G", "G": "T", "T": "A"}
	good := seqs["chr1"][999:1000]
	bad := next[seqs["chr1"][1999:2000]]
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t1000\t.\t" + good + "\t" + next[good] + "\t.\t.\tAF=1\n" +
		"chr1\t2000\t.\t" + bad + "\t" + next[bad] + "\t.\t.\tAF=1\n"
	if err := os.WriteFile(vcf, []byte(body), 0o600); err != nil {
		t.Fatal(err)
	}

	t.Run("mismatch warns", func(t *testing.T) {
		setGenerateFlags(t, map[string]string{"config": writePairedConfig(t, ""), "reference": ref, "stdout": "true", "variants": vcf, "seed": "1"})
		var out bytes.Buffer
		if err := runGenerate(generateCmd, &out); err != nil {
			t.Fatal(err)
		}
		if n := strings.Count(out.String(), "\n"); n != 2*2*200*4 {
			t.Errorf("%d lines", n)
		}
	})
	t.Run("strict-ref", func(t *testing.T) {
		setGenerateFlags(t, map[string]string{"config": writePairedConfig(t, ""), "reference": ref, "stdout": "true", "variants": vcf, "strict-ref": "true"})
		var out bytes.Buffer
		err := runGenerate(generateCmd, &out)
		if err == nil || !strings.Contains(err.Error(), "chr1:2000") {
			t.Fatalf("err = %v", err)
		}
		if out.Len() != 0 {
			t.Errorf("%d bytes written", out.Len())
		}
	})
	t.Run("hmm model", func(t *testing.T) {
		setGenerateFlags(t, map[string]string{"config": writePairedConfig(t, `,"hmm_model_path":"m.hmm.json"`), "reference": ref})
		if err := runGenerate(generateCmd, io.Discard); err == nil || !strings.Contains(err.Error(), "HMM model") {
			t.Fatalf("err = %v", err)
		}
	})
	t.Run("config reference_path", func(t *testing.T) {
		setGenerateFlags(t, map[string]string{"config": writePairedConfig(t, fmt.Sprintf(`,"reference_path":%q,"reference_region":"chr1"`, ref)), "stdout": "true", "total-reads": "3"})
		var out bytes.Buffer
		if err := runGenerate(generateCmd, &out); err != nil {
			t.Fatal(err)
		}
	})
	t.Run("gzip reference", func(t *testing.T) {
		gz := filepath.Join(t.TempDir(), "ref.fa.gz")
		if err := os.WriteFile(gz, []byte{0x1f, 0x8b, 8, 0}, 0o600); err != nil {
			t.Fatal(err)
		}
		out := filepath.Join(t.TempDir(), "out")
		setGenerateFlags(t, map[string]string{"config": writePairedConfig(t, ""), "reference": gz, "output": out})
		if err := runGenerate(generateCmd, io.Discard); err == nil || !strings.Contains(err.Error(), "gzip") {
			t.Fatalf("err = %v", err)
		}
		if _, err := os.Stat(out); !os.IsNotExist(err) {
			t.Errorf("output directory created before the reference error: %v", err)
		}
	})
}

// generate --variants in file mode (#32): with --reference the alleles are
// planted at their AF (AF 1: every fragment carrying the site has ALT), a
// variants table per batch lists them, the manifest records the VCF by
// name and sha256 (not its path), --strict-ref fails before writing, and
// variants_path in the config works like --variants.
func TestGenerate_FileModeVariants(t *testing.T) {
	ref, seqs := writeCmdReference(t)
	vcf := filepath.Join(t.TempDir(), "calls.vcf")
	next := map[string]string{"A": "C", "C": "G", "G": "T", "T": "A"}
	good := seqs["chr2"][2999:3000]
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr2\t3000\tsnv1\t" + good + "\t" + next[good] + "\t.\t.\tAF=1\n"
	if err := os.WriteFile(vcf, []byte(body), 0o600); err != nil {
		t.Fatal(err)
	}
	run := func(t *testing.T, flags map[string]string) string {
		out := filepath.Join(t.TempDir(), "out")
		flags["output"], flags["reference"], flags["region"], flags["illumina-errors"], flags["seed"] = out, ref, "chr2:2001-4000", "none", "3"
		setGenerateFlags(t, flags)
		if err := runGenerate(generateCmd, io.Discard); err != nil {
			t.Fatal(err)
		}
		return out
	}
	check := func(t *testing.T, out string) []byte {
		m, err := manifest.Read(filepath.Join(out, manifest.FileName))
		if err != nil {
			t.Fatal(err)
		}
		if err := m.VerifyDir(out); err != nil {
			t.Fatal(err)
		}
		n, sum, _ := manifest.HashFile(vcf)
		if m.Variants == nil || m.Variants.Name != "calls.vcf" || m.Variants.Bytes != n || m.Variants.SHA256 != sum || m.Config.VariantsPath != "" {
			t.Errorf("manifest variants %+v, config variants_path %q", m.Variants, m.Config.VariantsPath)
		}
		var all []byte
		alt := 0
		for _, o := range m.Outputs {
			if !strings.HasSuffix(o.Path, ".variants.tsv") {
				continue
			}
			data, err := os.ReadFile(filepath.Join(out, o.Path))
			if err != nil {
				t.Fatal(err)
			}
			all = append(all, data...)
			for _, line := range strings.Split(strings.TrimSpace(string(data)), "\n")[1:] {
				f := strings.Split(line, "\t")
				if f[4] != "snv1" || f[7] != "ALT" {
					t.Fatalf("variants line %q: want snv1 ALT", line)
				}
				alt++
			}
		}
		if alt < 20 {
			t.Errorf("only %d fragments carry the AF 1 SNV", alt)
		}
		return all
	}
	var byFlag, byConfig []byte
	t.Run("flag", func(t *testing.T) {
		byFlag = check(t, run(t, map[string]string{"config": writePairedConfig(t, ""), "variants": vcf}))
	})
	t.Run("config", func(t *testing.T) {
		byConfig = check(t, run(t, map[string]string{"config": writePairedConfig(t, fmt.Sprintf(`,"variants_path":%q`, vcf))}))
	})
	if len(byFlag) == 0 || !bytes.Equal(byFlag, byConfig) {
		t.Error("variants_path in the config and --variants plant different alleles")
	}
	t.Run("strict-ref", func(t *testing.T) {
		bad := filepath.Join(t.TempDir(), "bad.vcf")
		if err := os.WriteFile(bad, []byte(strings.Replace(body, "\t"+good+"\t", "\t"+next[good]+"\t", 1)), 0o600); err != nil {
			t.Fatal(err)
		}
		out := filepath.Join(t.TempDir(), "out")
		setGenerateFlags(t, map[string]string{"config": writePairedConfig(t, ""), "reference": ref, "variants": bad, "strict-ref": "true", "output": out})
		if err := runGenerate(generateCmd, io.Discard); err == nil || !strings.Contains(err.Error(), "chr2:3000") {
			t.Fatalf("err = %v", err)
		}
		if _, err := os.Stat(out); !os.IsNotExist(err) {
			t.Errorf("--strict-ref wrote output: %v", err)
		}
	})
}
