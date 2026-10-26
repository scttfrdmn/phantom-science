//go:build !release

package clinical

// Scenario reads sampled from a reference genome (--reference), and the
// platform read model of every scenario output (#24, #27).

import (
	"bufio"
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"math"
	mathrand "math/rand"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/illumina"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// fqRecord is one FASTQ record.
type fqRecord struct{ name, seq, qual string }

// parseRecords parses FASTQ data, failing on a malformed record.
func parseRecords(t *testing.T, data []byte) []fqRecord {
	t.Helper()
	sc := bufio.NewScanner(bytes.NewReader(data))
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	var lines []string
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	require.Zero(t, len(lines)%4, "FASTQ line count")
	var out []fqRecord
	for i := 0; i < len(lines); i += 4 {
		require.True(t, strings.HasPrefix(lines[i], "@"), lines[i])
		require.Equal(t, "+", lines[i+2])
		require.Len(t, lines[i+3], len(lines[i+1]))
		out = append(out, fqRecord{strings.SplitN(lines[i][1:], " ", 2)[0], lines[i+1], lines[i+3]})
	}
	return out
}

func readRecords(t *testing.T, path string) []fqRecord {
	t.Helper()
	data, err := os.ReadFile(path)
	require.NoError(t, err)
	return parseRecords(t, data)
}

// requireBinnedQualities checks every quality is one of the RTA3 bins.
func requireBinnedQualities(t *testing.T, recs []fqRecord) {
	t.Helper()
	ok := map[byte]bool{}
	for _, q := range illumina.Bins {
		ok['!'+q] = true
	}
	seen := map[byte]int{}
	for _, r := range recs {
		for i := 0; i < len(r.qual); i++ {
			require.True(t, ok[r.qual[i]], "quality %q is not a binned Illumina quality", r.qual[i])
			seen[r.qual[i]]++
		}
	}
	require.Greater(t, seen['F'], 0, "Q37 is the most common bin")
}

func rc(s string) string {
	b := []byte(s)
	for i, j := 0, len(b)-1; i <= j; i, j = i+1, j-1 {
		b[i], b[j] = comp(b[j]), comp(b[i])
	}
	return string(b)
}

func comp(c byte) byte {
	switch c {
	case 'A':
		return 'T'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'T':
		return 'A'
	}
	return c
}

// seqIndex locates reads in a few sequences by their first k-mer.
type seqIndex struct {
	seqs []string
	k    int
	at   map[string][][2]int // k-mer -> (sequence, offset)
}

func newSeqIndex(k int, seqs ...[]byte) *seqIndex {
	ix := &seqIndex{k: k, at: map[string][][2]int{}}
	for s, b := range seqs {
		ix.seqs = append(ix.seqs, string(b))
		for i := 0; i+k <= len(b); i++ {
			key := string(b[i : i+k])
			ix.at[key] = append(ix.at[key], [2]int{s, i})
		}
	}
	return ix
}

// locate returns where read occurs exactly, forward (minus false) or as
// its reverse complement: the sequence index and the 0-based forward
// offset of its leftmost base.
func (ix *seqIndex) locate(read string) (seq, off int, minus, ok bool) {
	for _, strand := range []bool{false, true} {
		r := read
		if strand {
			r = rc(read)
		}
		if len(r) < ix.k {
			continue
		}
		for _, h := range ix.at[r[:ix.k]] {
			if strings.HasPrefix(ix.seqs[h[0]][h[1]:], r) {
				return h[0], h[1], strand, true
			}
		}
	}
	return 0, 0, false, false
}

// countCovering counts the reads containing each context, on either strand.
func countCovering(recs []fqRecord, contexts ...string) []int {
	n := make([]int, len(contexts))
	for _, r := range recs {
		fwd, rev := r.seq, rc(r.seq)
		for i, c := range contexts {
			if strings.Contains(fwd, c) || strings.Contains(rev, c) {
				n[i]++
			}
		}
	}
	return n
}

// requireWithin4Sigma checks that k of n draws is within 4σ of n·p.
func requireWithin4Sigma(t *testing.T, k, n int, p float64, what string) {
	t.Helper()
	require.Greater(t, n, 30, "%s: too few reads cover the site", what)
	sd := math.Sqrt(float64(n) * p * (1 - p))
	require.LessOrEqual(t, math.Abs(float64(k)-float64(n)*p), 4*sd,
		"%s: %d ALT of %d, expected %.1f ± 4×%.1f", what, k, n, float64(n)*p, sd)
	t.Logf("%s: %d/%d ALT (VAF %.3f, expected %.2f)", what, k, n, float64(k)/float64(n), p)
}

// writeRandomFASTA writes a one-contig reference "chr1" of n random bases,
// 70 bases per line, with bases 2000-2499 soft-masked, and returns its path
// and the upper-case sequence.
func writeRandomFASTA(t *testing.T, n int, seed int64) (string, []byte) {
	t.Helper()
	// #nosec G404 -- test fixture
	rng := mathrand.New(mathrand.NewSource(seed))
	seq := make([]byte, n)
	for i := range seq {
		seq[i] = "ACGT"[rng.Intn(4)]
	}
	file := append([]byte(nil), seq...)
	copy(file[2000:2500], bytes.ToLower(seq[2000:2500]))
	var buf bytes.Buffer
	buf.WriteString(">chr1 fixture\n")
	for a := 0; a < n; a += 70 {
		buf.Write(file[a:min(a+70, n)])
		buf.WriteByte('\n')
	}
	path := filepath.Join(t.TempDir(), "ref.fa")
	require.NoError(t, os.WriteFile(path, buf.Bytes(), 0o600))
	return path, seq
}

// refScenario is a scenario run against a reference fixture.
type refScenario struct {
	spec    PlantedVariantSpec
	paired  bool
	trio    bool
	errors  string // illumina errors; "" is the default (quality)
	seed    uint64
	cov     int
	readLen int
}

// run runs the scenario with reference ref (or synthetic bases for "")
// and returns the manifest and the output paths by member ("proband",
// "proband_r2", ...).
func (rs refScenario) run(t *testing.T, ref string) (*GroundTruthManifest, map[string]string) {
	t.Helper()
	dir := t.TempDir()
	out := map[string]string{"proband": filepath.Join(dir, "proband_R1.fastq")}
	sc := &ClinicalScenario{
		Phenotype: PhenotypeSpec{PlantedVariant: &rs.spec},
		Sequencer: SequencerSpec{Platform: "NovaSeq6000", Coverage: rs.cov, ReadLength: rs.readLen, PairedEnd: rs.paired},
		Output:    OutputSpec{ProbandFASTQ: out["proband"], GroundTruth: filepath.Join(dir, "truth.json")},
	}
	if rs.trio {
		sc.Trio = &TrioSpec{Generate: true, DeNovo: true}
		out["mother"] = filepath.Join(dir, "mother_R1.fastq")
		out["father"] = filepath.Join(dir, "father_R1.fastq")
		sc.Output.MotherFASTQ, sc.Output.FatherFASTQ = out["mother"], out["father"]
	}
	if rs.paired {
		out["proband_r2"] = filepath.Join(dir, "proband_R2.fastq")
		sc.Output.ProbandFASTQR2 = out["proband_r2"]
		if rs.trio {
			out["mother_r2"] = filepath.Join(dir, "mother_R2.fastq")
			out["father_r2"] = filepath.Join(dir, "father_R2.fastq")
			sc.Output.MotherFASTQR2, sc.Output.FatherFASTQR2 = out["mother_r2"], out["father_r2"]
		}
	}
	cfg := rs.config()
	g := NewScenarioGenerator(&cfg)
	if ref != "" {
		g.WithReference(ref)
	}
	m, err := g.Run(context.Background(), sc)
	require.NoError(t, err)
	require.NoError(t, WriteManifest(sc.Output.GroundTruth, m))
	return m, out
}

func (rs refScenario) config() types.SequencerConfig {
	return types.SequencerConfig{SampleName: "ref", Lane: 1, ReadLength: rs.readLen, Seed: rs.seed,
		Population: "CEU", Illumina: types.IlluminaSpec{Errors: rs.errors}}
}

// applyVariant returns seq (0-based coordinates) with the 1-based variant
// applied.
func applyVariant(seq []byte, v PlantedVariantSpec) []byte {
	p := int(v.Pos - 1)
	out := append([]byte(nil), seq[:p]...)
	out = append(out, v.Alt...)
	return append(out, seq[p+len(v.Ref):]...)
}

// Proband reads with a reference: with sequencing errors off every read is
// an exact substring (or reverse complement) of the reference, or of the
// reference with the planted variant; the planted ALT appears at the
// expected VAF (within 4σ of the binomial over the reads covering it); the
// qualities are binned; the output is a function of the seed.
func TestScenarioReference_ProbandReadsFromReference(t *testing.T) {
	ref, seq := writeRandomFASTA(t, 12000, 1)
	snvRef := string(seq[6000])
	snvAlt := string("ACGT"[(strings.IndexByte("ACGT", seq[6000])+1)%4])
	cases := []struct {
		name string
		spec PlantedVariantSpec
		vaf  float64
	}{
		{"het SNV", PlantedVariantSpec{Chrom: "chr1", Pos: 6001, Ref: snvRef, Alt: snvAlt, Zygosity: "heterozygous"}, 0.5},
		{"hom SNV", PlantedVariantSpec{Chrom: "1", Pos: 6001, Ref: snvRef, Alt: snvAlt, Zygosity: "homozygous"}, 1},
		{"het deletion", PlantedVariantSpec{Chrom: "chr1", Pos: 7001, Ref: string(seq[7000:7004]), Alt: string(seq[7000]), Zygosity: "heterozygous"}, 0.5},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			rs := refScenario{spec: c.spec, errors: types.IlluminaErrorsNone, seed: 11, cov: 100, readLen: 100}
			m, out := rs.run(t, ref)
			assert.Equal(t, basesReference, m.Sequencer.Bases)
			want, err := manifest.HashInput(ref)
			require.NoError(t, err)
			assert.Equal(t, want, m.Reference)
			assert.True(t, strings.HasSuffix(m.Sequencer.RegionFocused, "1:0-12000"), "window clipped to the contig: %s", m.Sequencer.RegionFocused)
			recs := readRecords(t, out["proband"])
			require.Equal(t, 100*12000/100, len(recs))
			assert.Equal(t, len(recs), m.ComputeMetrics.ReadsGenerated)
			requireBinnedQualities(t, recs)

			alt := applyVariant(seq, c.spec)
			ix := newSeqIndex(24, seq, alt)
			for _, r := range recs {
				require.Len(t, r.seq, 100)
				_, _, _, ok := ix.locate(r.seq)
				require.True(t, ok, "read %s is not in the reference or its ALT haplotype: %s", r.name, r.seq)
			}
			p := int(c.spec.Pos - 1)
			refCtx := string(seq[p-8 : p+len(c.spec.Ref)+8])
			altCtx := string(alt[p-8 : p+len(c.spec.Alt)+8])
			n := countCovering(recs, altCtx, refCtx)
			requireWithin4Sigma(t, n[0], n[0]+n[1], c.vaf, c.name)
			if c.vaf == 1 {
				assert.Zero(t, n[1], "homozygous: no REF reads")
			}

			// Determinism.
			_, out2 := rs.run(t, ref)
			a, err := os.ReadFile(out["proband"])
			require.NoError(t, err)
			b, err := os.ReadFile(out2["proband"])
			require.NoError(t, err)
			assert.True(t, bytes.Equal(a, b), "same seed, same FASTQ")
			rs.seed++
			_, out3 := rs.run(t, ref)
			d, err := os.ReadFile(out3["proband"])
			require.NoError(t, err)
			assert.False(t, bytes.Equal(a, d), "another seed, other reads")
		})
	}
}

// With the default error model the reads carry quality-matched
// substitutions: most reads still match the reference exactly, some do not.
func TestScenarioReference_DefaultErrorModel(t *testing.T) {
	ref, seq := writeRandomFASTA(t, 12000, 2)
	spec := PlantedVariantSpec{Chrom: "chr1", Pos: 6001, Ref: string(seq[6000]), Alt: string(comp(seq[6000])), Zygosity: "heterozygous"}
	_, out := refScenario{spec: spec, seed: 5, cov: 30, readLen: 150}.run(t, ref)
	recs := readRecords(t, out["proband"])
	requireBinnedQualities(t, recs)
	ix := newSeqIndex(24, seq, applyVariant(seq, spec))
	exact := 0
	for _, r := range recs {
		if _, _, _, ok := ix.locate(r.seq); ok {
			exact++
		}
	}
	frac := float64(exact) / float64(len(recs))
	t.Logf("%d of %d reads (%.1f%%) match the reference exactly", exact, len(recs), 100*frac)
	assert.Greater(t, frac, 0.3)
	assert.Less(t, frac, 0.999, "sequencing errors are applied")
}

// Paired-end scenario reads are the two ends of one reference fragment:
// R1 forward and R2 reverse-complemented (or the reverse), the insert
// between the read length and 1 kb, written to separate R1 and R2 files with
// matching names.
func TestScenarioReference_PairedFragments(t *testing.T) {
	ref, seq := writeRandomFASTA(t, 12000, 3)
	spec := PlantedVariantSpec{Chrom: "chr1", Pos: 6001, Ref: string(seq[6000]), Alt: string(comp(seq[6000])), Zygosity: "heterozygous"}
	m, out := refScenario{spec: spec, paired: true, errors: types.IlluminaErrorsNone, seed: 7, cov: 30, readLen: 100}.run(t, ref)
	assert.True(t, m.Sequencer.PairedEnd)
	r1, r2 := readRecords(t, out["proband"]), readRecords(t, out["proband_r2"])
	require.Equal(t, 30*12000/200, len(r1), "pairs for 30x")
	require.Equal(t, len(r1), len(r2))
	assert.Equal(t, 2*len(r1), m.ComputeMetrics.ReadsGenerated)
	requireBinnedQualities(t, append(r1, r2...))
	ix := newSeqIndex(24, seq, applyVariant(seq, spec))
	var inserts []float64
	for i := range r1 {
		require.Equal(t, r1[i].name, r2[i].name, "mates share a name")
		// An SNV keeps coordinates, so offsets in the reference and in its
		// ALT haplotype agree.
		_, o1, m1, ok1 := ix.locate(r1[i].seq)
		_, o2, m2, ok2 := ix.locate(r2[i].seq)
		require.True(t, ok1 && ok2, "pair %d not in the reference", i)
		require.NotEqual(t, m1, m2, "mates are on opposite strands")
		insert := o2 + 100 - o1
		if m1 {
			insert = o1 + 100 - o2
		}
		require.GreaterOrEqual(t, insert, 100)
		require.LessOrEqual(t, insert, 1000)
		inserts = append(inserts, float64(insert))
	}
	var mean float64
	for _, x := range inserts {
		mean += x
	}
	mean /= float64(len(inserts))
	assert.InDelta(t, types.DefaultInsertSizeMean, mean, 10, "mean insert")
}

// Trio scenarios, with and without a reference: each member's haplotypes
// differ from the window bases only at the simulated background sites (and
// the de novo site on the proband's first haplotype), with the alleles the
// simulation gave each haplotype; the proband's alleles come from its
// parents; every read (errors off) is an exact substring of one of the
// member's haplotypes; the de novo ALT is at VAF 1/2 in the proband and
// absent from the parents; qualities are binned; and the FASTQs are a
// function of the seed.
func TestScenarioReference_Trio(t *testing.T) {
	ref, seq := writeRandomFASTA(t, 20000, 4)
	spec := PlantedVariantSpec{Chrom: "chr1", Pos: 10001, Ref: string(seq[10000]), Alt: string(comp(seq[10000])), Zygosity: "heterozygous"}
	for _, withRef := range []bool{true, false} {
		t.Run(fmt.Sprintf("reference=%v", withRef), func(t *testing.T) {
			path := ""
			if withRef {
				path = ref
			}
			rs := refScenario{spec: spec, trio: true, errors: types.IlluminaErrorsNone, seed: 2026, cov: 60, readLen: 100}
			m, out := rs.run(t, path)
			require.NotNil(t, m.Trio)
			want := basesSynthetic
			if withRef {
				want = basesReference
			}
			assert.Equal(t, want, m.Sequencer.Bases)

			// The plan the run used (it is a function of the seed).
			cfg := rs.config()
			g := NewScenarioGenerator(&cfg)
			region := &genomics.GenomicRegion{Chromosome: "chr1", Start: 0, End: 35000}
			if withRef {
				g.WithReference(ref)
				region.End = 20000
			}
			dn := &genomics.VariantRecord{Chromosome: "chr1", Position: spec.Pos, Ref: spec.Ref, Alt: spec.Alt, Freq: 0.5}
			tp, err := g.planTrio(cfg, region, "chr1", dn)
			require.NoError(t, err)
			if withRef {
				require.Equal(t, string(seq), string(tp.base), "the window is the reference")
			}
			require.Equal(t, m.Trio.BackgroundVariants, tp.sim.NumVariants)

			haps := map[string][2][]byte{}
			for _, member := range []string{"proband", "mother", "father"} {
				h, err := tp.memberHaplotypes(member)
				require.NoError(t, err)
				haps[member] = h
			}
			bits := map[string][2][]int{
				"proband": {tp.trio.ProbandH1, tp.trio.ProbandH2},
				"mother":  {tp.trio.MotherH1, tp.trio.MotherH2},
				"father":  {tp.trio.FatherH1, tp.trio.FatherH2},
			}
			sites := map[int]int{} // window offset -> variant index
			for i, v := range tp.trio.Variants {
				require.Len(t, v.Alleles, 2, "every site on a fixture base is used")
				require.Equal(t, tp.base[v.Position], v.Alleles[0][0], "site %d: REF is the window base", v.Position)
				require.NotEqual(t, v.Alleles[0], v.Alleles[1], "site %d: ALT differs from REF", v.Position)
				simAlleles := simulatedAlleles(t, tp, v.Position)
				if !withRef {
					require.Equal(t, simAlleles, v.Alleles, "synthetic: the simulated alleles")
				} else {
					require.Contains(t, simAlleles, v.Alleles[1], "reference: ALT is a simulated base")
				}
				sites[v.Position] = i
			}
			require.Greater(t, len(sites), 5, "background sites in the window")
			dn0 := int(spec.Pos - 1)
			for member, hs := range haps {
				for h, hap := range hs {
					require.Len(t, hap, len(tp.base))
					for p := range hap {
						i, isSite := sites[p]
						switch {
						case p == dn0:
							wantBase := spec.Ref[0]
							if member == "proband" && h == 0 {
								wantBase = spec.Alt[0]
							}
							require.Equal(t, wantBase, hap[p], "%s h%d de novo site", member, h+1)
						case isSite:
							require.Equal(t, tp.trio.Variants[i].Alleles[bits[member][h][i]][0], hap[p], "%s h%d site %d", member, h+1, p)
						default:
							require.Equal(t, tp.base[p], hap[p], "%s h%d: base %d differs from the window", member, h+1, p)
						}
					}
				}
			}
			// Mendelian: proband H1 maternal, H2 paternal, at every site.
			for p := range sites {
				assert.Contains(t, []byte{haps["mother"][0][p], haps["mother"][1][p]}, haps["proband"][0][p], "site %d maternal", p)
				assert.Contains(t, []byte{haps["father"][0][p], haps["father"][1][p]}, haps["proband"][1][p], "site %d paternal", p)
			}

			altCtx := string(haps["proband"][0][dn0-8 : dn0+9])
			for _, member := range []string{"proband", "mother", "father"} {
				recs := readRecords(t, out[member])
				require.Equal(t, 60*len(tp.base)/100, len(recs), member)
				requireBinnedQualities(t, recs)
				ix := newSeqIndex(24, haps[member][0], haps[member][1])
				for _, r := range recs {
					_, _, _, ok := ix.locate(r.seq)
					require.True(t, ok, "%s read %s is not on its haplotypes", member, r.name)
				}
				refCtx := string(haps[member][1][dn0-8 : dn0+9])
				n := countCovering(recs, altCtx, refCtx)
				if member == "proband" {
					requireWithin4Sigma(t, n[0], n[0]+n[1], 0.5, "proband de novo")
				} else {
					assert.Zero(t, n[0], "%s has no de novo ALT", member)
				}
			}

			_, out2 := rs.run(t, path)
			for _, member := range []string{"proband", "mother", "father"} {
				a, err := os.ReadFile(out[member])
				require.NoError(t, err)
				b, err := os.ReadFile(out2[member])
				require.NoError(t, err)
				assert.True(t, bytes.Equal(a, b), "%s: same seed, same FASTQ", member)
			}
			a, _ := os.ReadFile(out["mother"])
			b, _ := os.ReadFile(out["father"])
			assert.False(t, bytes.Equal(a, b), "parents draw independent reads")
		})
	}
}

// simulatedAlleles returns the coalescent's alleles at pos.
func simulatedAlleles(t *testing.T, tp *trioPlan, pos int) []string {
	t.Helper()
	for _, v := range tp.sim.Variants {
		if v.Position == pos {
			return v.Alleles
		}
	}
	t.Fatalf("no simulated site at %d", pos)
	return nil
}

// The manifest records the read model and the reference by name, size and
// SHA-256 (#38). It holds no local path and no wall-clock value: two runs
// of the same scenario and seed, writing to different directories, write
// the same bytes.
func TestScenarioReference_ManifestJSON(t *testing.T) {
	ref, seq := writeRandomFASTA(t, 12000, 5)
	spec := PlantedVariantSpec{Chrom: "chr1", Pos: 6001, Ref: string(seq[6000]), Alt: string(comp(seq[6000])), Zygosity: "heterozygous"}
	rs := refScenario{spec: spec, paired: true, seed: 3, cov: 2, readLen: 100}
	var runs [2][]byte
	var dirs [2]string
	for i := range runs {
		_, out := rs.run(t, ref)
		dirs[i] = filepath.Dir(out["proband"])
		data, err := os.ReadFile(filepath.Join(dirs[i], "truth.json"))
		require.NoError(t, err)
		runs[i] = data
	}
	assert.Equal(t, string(runs[0]), string(runs[1]), "the same scenario and seed wrote different manifests")
	for _, p := range []string{ref, filepath.Dir(ref), dirs[0], os.TempDir()} {
		assert.NotContains(t, string(runs[0]), p, "the manifest holds a local path")
	}
	var m map[string]any
	require.NoError(t, json.Unmarshal(runs[0], &m))
	sq := m["sequencer"].(map[string]any)
	assert.Equal(t, true, sq["paired_end"])
	assert.Equal(t, "reference", sq["bases"])
	assert.NotContains(t, m, "generated_at")
	assert.NotContains(t, m["compute_metrics"], "wall_clock_ns")
	n, sum, err := manifest.HashFile(ref)
	require.NoError(t, err)
	assert.Equal(t, map[string]any{"name": "ref.fa", "bytes": float64(n), "sha256": sum}, m["reference"])
}
