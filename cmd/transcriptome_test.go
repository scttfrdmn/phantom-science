//go:build !release

package cmd

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

const txFixtures = "../internal/genomics/transcriptome/testdata"

var transcriptomeFlagNames = []string{
	"transcripts", "genome", "gtf", "expression-profile", "reads", "out", "tier", "seed",
	"novel-rate", "compress", "kit", "full-length-fraction", "zipf-alpha", "expression-sigma", "sample-id",
}

// runTranscriptomeCLI runs "phantom-science transcriptome args..." through
// the root command and resets the subcommand's flags afterwards.
func runTranscriptomeCLI(t *testing.T, args ...string) (stdout, stderr string, err error) {
	t.Helper()
	var out, errb bytes.Buffer
	rootCmd.SetOut(&out)
	rootCmd.SetErr(&errb)
	rootCmd.SetArgs(append([]string{"transcriptome"}, args...))
	defer func() {
		rootCmd.SetArgs([]string{})
		rootCmd.SetOut(nil)
		rootCmd.SetErr(nil)
		for _, name := range transcriptomeFlagNames {
			f := transcriptomeCmd.Flags().Lookup(name)
			_ = f.Value.Set(f.DefValue)
			f.Changed = false
		}
	}()
	err = rootCmd.Execute()
	return out.String(), errb.String(), err
}

func TestTranscriptomeCommand_Flags(t *testing.T) {
	for _, name := range transcriptomeFlagNames {
		if transcriptomeCmd.Flags().Lookup(name) == nil {
			t.Errorf("missing --%s", name)
		}
	}
}

func TestTranscriptomeCommand_GenomeGTF(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "run")
	args := []string{
		"--genome", filepath.Join(txFixtures, "tiny_genome.fna"), "--gtf", filepath.Join(txFixtures, "tiny.gtf"),
		"--reads", "200", "--out", dir, "--seed", "5", "--novel-rate", "0.3", "--compress", "gzip",
	}
	stdout, stderr, err := runTranscriptomeCLI(t, args...)
	if err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	if !strings.Contains(stdout, "Wrote 200 cDNA reads") || !strings.Contains(stderr, "1 annotated transcripts are on chromosomes missing") {
		t.Errorf("stdout %q stderr %q", stdout, stderr)
	}
	for _, f := range []string{"reads.fastq.gz", "sequencing_summary.txt", "truth.tsv", "counts.tsv", "novel_isoforms.tsv", "novel_isoforms.gtf"} {
		if st, err := os.Stat(filepath.Join(dir, f)); err != nil || st.Size() == 0 {
			t.Errorf("%s missing or empty", f)
		}
	}
	truth, _ := os.ReadFile(filepath.Join(dir, "truth.tsv"))
	if n := strings.Count(string(truth), "\n"); n != 201 {
		t.Errorf("truth.tsv has %d lines, want 201", n)
	}
	counts, _ := os.ReadFile(filepath.Join(dir, "counts.tsv"))
	total := 0
	for _, l := range strings.Split(strings.TrimSpace(string(counts)), "\n")[1:] {
		f := strings.Split(l, "\t")
		n, _ := strconv.Atoi(f[len(f)-1])
		total += n
	}
	if total != 200 {
		t.Errorf("counts.tsv reads sum to %d", total)
	}

	// Same seed: byte-identical output.
	dir2 := filepath.Join(t.TempDir(), "run")
	args[len(args)-7] = dir2
	if _, stderr, err := runTranscriptomeCLI(t, args...); err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	for _, f := range []string{"reads.fastq.gz", "truth.tsv", "sequencing_summary.txt"} {
		a, _ := os.ReadFile(filepath.Join(dir, f))
		b, _ := os.ReadFile(filepath.Join(dir2, f))
		if !bytes.Equal(a, b) {
			t.Errorf("%s differs between runs with --seed 5", f)
		}
	}
}

func TestTranscriptomeCommand_TranscriptsAndProfile(t *testing.T) {
	tmp := t.TempDir()
	fa := filepath.Join(tmp, "tx.fa")
	os.WriteFile(fa, []byte(">t1|g1|\n"+strings.Repeat("ACGTTGCA", 40)+"\n>t2.3|g2|\n"+strings.Repeat("GGCATTCA", 30)+"\n"), 0o600)
	prof := filepath.Join(tmp, "tpm.tsv")
	os.WriteFile(prof, []byte("transcript_id\ttpm\nt1\t0\nt2\t5\nzz\t1\n"), 0o600)
	out := filepath.Join(tmp, "o")
	_, stderr, err := runTranscriptomeCLI(t, "--transcripts", fa, "--expression-profile", prof,
		"--reads", "50", "--out", out, "--tier", "none")
	if err != nil {
		t.Fatalf("%v\n%s", err, stderr)
	}
	if !strings.Contains(stderr, "[info] seed ") || !strings.Contains(stderr, "2 matched, 1 not in the transcript table") {
		t.Errorf("stderr %q lacks the drawn seed / profile summary", stderr)
	}
	truth, _ := os.ReadFile(filepath.Join(out, "truth.tsv"))
	for _, l := range strings.Split(strings.TrimSpace(string(truth)), "\n")[1:] {
		if f := strings.Split(l, "\t"); f[1] != "t2.3" {
			t.Fatalf("read from %s, but only t2 (version-stripped match) is expressed", f[1])
		}
	}
	if _, err := os.Stat(filepath.Join(out, "novel_isoforms.gtf")); err == nil {
		t.Error("novel GTF written without exon coordinates")
	}
}

func TestTranscriptomeCommand_Errors(t *testing.T) {
	fa := filepath.Join(txFixtures, "tiny_genome.fna")
	gtf := filepath.Join(txFixtures, "tiny.gtf")
	out := t.TempDir()
	for _, c := range []struct {
		want string
		args []string
	}{
		{"--out is required", []string{"--transcripts", fa, "--reads", "1"}},
		{"--reads must be positive", []string{"--transcripts", fa, "--out", out}},
		{"not both", []string{"--transcripts", fa, "--genome", fa, "--gtf", gtf, "--reads", "1", "--out", out}},
		{"--genome requires --gtf", []string{"--genome", fa, "--reads", "1", "--out", out}},
		{"is required", []string{"--reads", "1", "--out", out}},
		{"unknown --compress", []string{"--transcripts", fa, "--reads", "1", "--out", out, "--compress", "bz2"}},
		{"unknown basecaller tier", []string{"--transcripts", fa, "--reads", "1", "--out", out, "--tier", "ultra"}},
		{"unknown cDNA kit", []string{"--transcripts", fa, "--reads", "1", "--out", out, "--kit", "PCS999"}},
		{"novel-isoform rate", []string{"--transcripts", fa, "--reads", "1", "--out", out, "--novel-rate", "2"}},
	} {
		if _, _, err := runTranscriptomeCLI(t, c.args...); err == nil || !strings.Contains(err.Error(), c.want) {
			t.Errorf("%v: err = %v, want %q", c.args, err, c.want)
		}
	}
}
