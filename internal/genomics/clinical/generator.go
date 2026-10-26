//go:build !release

package clinical

// How scenario reads are made.
//
// Every scenario FASTQ is written by the engine's own read pipeline
// (Engine.GenerateToWriter), so scenario reads have the platform quality and
// error model of `generate`: for Illumina, binned qualities, quality-matched
// substitution errors (illumina.errors), and with sequencer.paired_end the
// two ends of one fragment (insert sizes from illumina.insert_size).
//
// Proband-only scenarios: with a reference, the engine samples fragments
// from the reference over the scenario window (reference_region) and plants
// the variant at its coordinate with probability ExpectedVAF per fragment;
// without one, reads are synthetic bases positioned in the window
// (GenerateOpts.RefRegion), with the variant planted the same way.
//
// Trio scenarios: the window's bases (the reference, or one synthetic
// sequence drawn from the engine) carry the simulated background variants;
// each member's two haplotypes are built from them (the proband's first
// haplotype also carries the de novo variant) and written as a two-contig
// FASTA, and the engine samples that member's reads from it. Each read
// therefore comes from one haplotype, with every site on it phased.

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"io"
	mathrand "math/rand"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	runmanifest "github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"github.com/scttfrdmn/phantom-science/internal/population"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

const windowHalfSize = 25000 // 50 kb window around variant (±25 kb)

// trioSimSeedSalt decorrelates the haplotype simulation's RNG from the
// engines' read streams, which use the same run seed.
const trioSimSeedSalt = 0x5452494f48415053 // "TRIOHAPS"

// backgroundStream is the engine stream the synthetic trio window is drawn
// from (without a reference).
const backgroundStream = 4

// Read-base sources recorded in the manifest.
const (
	basesReference = "reference"
	basesSynthetic = "synthetic"
)

// ScenarioGenerator orchestrates haplotype-aware clinical scenario generation.
type ScenarioGenerator struct {
	config        *types.SequencerConfig
	referencePath string
	// observe, when set, receives the bases of every read instead of the
	// reads being written (Regenerate).
	observe func(member string, seq []byte) error
	// rec describes the last run, for its run manifest.
	rec *runRecord
}

// NewScenarioGenerator creates a ScenarioGenerator. A reference_path in
// config is used as the scenario reference unless WithReference sets one.
func NewScenarioGenerator(config *types.SequencerConfig) *ScenarioGenerator {
	return &ScenarioGenerator{config: config, referencePath: config.ReferencePath}
}

// WithReference makes Run take read bases from the reference FASTA at path
// (uncompressed; a samtools .fai next to it is used when present). The
// planted variant's REF allele is checked against it first (see
// CheckPlantedRef), and a mismatch fails the run before any reads are
// written. Reads are sampled from the scenario window of the reference;
// trio background variants are placed on the reference bases.
func (g *ScenarioGenerator) WithReference(path string) *ScenarioGenerator {
	g.referencePath = path
	return g
}

// output is one member's FASTQ destination: r1, and r2 for separate mate
// files ("" writes pairs interleaved in r1).
type output struct{ r1, r2 string }

// Run executes a clinical scenario and returns the ground truth manifest.
func (g *ScenarioGenerator) Run(ctx context.Context, scenario *ClinicalScenario) (*GroundTruthManifest, error) {
	start := time.Now()

	g.rec = nil
	var refInput *runmanifest.Input
	if g.referencePath != "" && g.observe == nil { // regeneration has checked the file
		var err error
		if refInput, err = runmanifest.HashInput(g.referencePath); err != nil {
			return nil, fmt.Errorf("reference: %w", err)
		}
	}

	// 1. Resolve variant.
	var spec *PlantedVariantSpec
	var gene string
	var entry *CatalogEntry
	if scenario.Phenotype.PlantedVariant != nil {
		spec = scenario.Phenotype.PlantedVariant
		gene = scenario.Phenotype.ResolvedGene
	} else {
		var err error
		spec, entry, err = ResolveHPO(scenario.Phenotype.HPOTerms)
		if err != nil {
			return nil, fmt.Errorf("resolving HPO terms: %w", err)
		}
		gene = entry.Gene
	}
	if g.referencePath != "" {
		if err := CheckPlantedRef(g.referencePath, spec); err != nil {
			return nil, err
		}
	}

	// 2. The ±25 kb window around the variant (inside the contig, with a
	// reference).
	region, refContig, err := g.scenarioRegion(spec)
	if err != nil {
		return nil, err
	}

	readLength := g.config.ReadLength
	if scenario.Sequencer.ReadLength > 0 {
		readLength = scenario.Sequencer.ReadLength
	}
	if readLength <= 0 {
		return nil, fmt.Errorf("scenario read length must be positive, got %d", readLength)
	}
	coverage := scenario.Sequencer.Coverage
	if coverage <= 0 {
		coverage = 30
	}
	paired := scenario.Sequencer.PairedEnd || g.config.PairedEnd
	mates := 1
	if paired {
		mates = 2
	}

	// Reads (pairs) for the requested coverage of the window.
	regionLen := region.End - region.Start
	units := int(int64(coverage) * regionLen / int64(readLength*mates))
	if units < 1 {
		units = 1
	}

	// 3. Build VariantRecord.
	vaf := spec.ExpectedVAF()
	varRecord := genomics.VariantRecord{
		Ref:        spec.Ref,
		Alt:        spec.Alt,
		Freq:       vaf,
		Chromosome: spec.Chrom,
		Position:   spec.Pos,
	}

	// Clone config for generation. The scenario decides where bases come
	// from, so a reference in the config is not passed through as is.
	cfg := *g.config
	cfg.ReadLength = readLength
	cfg.PairedEnd = paired
	cfg.ReadsPerBatch = units
	cfg.TotalBatches = 1
	cfg.ReferencePath, cfg.ReferenceRegion, cfg.StrictRef = "", "", false

	// Resolve the seed once so the haplotype simulation and every trio
	// member's engine derive from it; the scenario is reproducible with it.
	// With a keyed-seed fingerprint it is the seed derived from the key,
	// with the configured (or a drawn) seed as the salt: everything below
	// follows from it, so the scenario is the unfingerprinted scenario with
	// that seed (see fingerprint.go).
	rs, err := resolveScenarioSeed(cfg)
	if err != nil {
		return nil, err
	}
	switch {
	case g.observe != nil:
	case rs.Drawn && rs.Fingerprint != nil:
		fmt.Fprintf(os.Stderr, "   [info] fingerprint salt %d (keyed seed; set seed: %d in the sequencer config, with the same key, customer and dataset, to reproduce this scenario)\n",
			rs.Fingerprint.Salt, rs.Fingerprint.Salt)
	case rs.Drawn:
		fmt.Fprintf(os.Stderr, "   [info] seed %d (set seed: %d in the sequencer config to reproduce this scenario)\n", rs.Seed, rs.Seed)
	}
	cfg.Seed = rs.Seed
	if rs.Fingerprint != nil {
		// The member engines use the derived seed as it is.
		cfg.Fingerprinting, cfg.Enterprise.Fingerprinting = types.FingerprintConfig{}, types.FingerprintConfig{}
	}
	if g.observe != nil {
		cfg.Output.Compression = "" // reads are parsed, not written
	}
	rec := &runRecord{seed: rs, reference: refInput}

	outs := map[string]output{
		"proband": {scenario.Output.ProbandFASTQ, scenario.Output.ProbandFASTQR2},
		"mother":  {scenario.Output.MotherFASTQ, scenario.Output.MotherFASTQR2},
		"father":  {scenario.Output.FatherFASTQ, scenario.Output.FatherFASTQR2},
	}

	var readsGenerated int
	var simResult *population.HaplotypeSimulationResult

	if scenario.Trio != nil && scenario.Trio.Generate && scenario.Trio.DeNovo {
		// 4. Trio path: haplotype-aware generation.
		tp, err := g.planTrio(cfg, region, refContig, &varRecord)
		if err != nil {
			return nil, err
		}
		simResult = tp.sim
		for _, member := range []string{"proband", "mother", "father"} {
			if outs[member].r1 == "" {
				continue
			}
			n, err := g.writeMemberReads(cfg, tp, member, outs[member], rec)
			if err != nil {
				return nil, fmt.Errorf("writing %s FASTQ: %w", member, err)
			}
			readsGenerated += n
		}
	} else {
		// 5. Proband only: the engine plants the variant by coordinate.
		opts := genomics.GenerateOpts{Variants: []genomics.VariantRecord{varRecord}}
		if g.referencePath != "" {
			cfg.ReferencePath = g.referencePath
			cfg.ReferenceRegion = fmt.Sprintf("%s:%d-%d", refContig, region.Start+1, region.End)
			cfg.StrictRef = true
		} else {
			opts.RefRegion = region
		}
		n, err := g.writeReads(cfg, opts, "proband", outs["proband"], rec)
		if err != nil {
			return nil, fmt.Errorf("generating proband reads: %w", err)
		}
		readsGenerated = n
	}

	if g.observe == nil {
		if err := rec.finish(g, scenario, cfg); err != nil {
			return nil, err
		}
		g.rec = rec
	}

	// 6. Build manifest.
	planted := PlantedVariant{
		Gene:        gene,
		Chrom:       spec.Chrom,
		Pos:         spec.Pos,
		Ref:         spec.Ref,
		Alt:         spec.Alt,
		Zygosity:    spec.Zygosity,
		ExpectedVAF: vaf,
		Consequence: spec.Consequence,
		ACMGClass:   spec.ACMGClass,
		HPOTerms:    scenario.Phenotype.HPOTerms,
	}
	if entry != nil {
		planted.Build = entry.Build
		planted.ClinVarID = entry.ClinVarID
		planted.HGVS = entry.HGVS
	}
	bases := basesSynthetic
	if g.referencePath != "" {
		bases = basesReference
	}
	manifest := &GroundTruthManifest{
		ProbandID:       g.config.SampleName,
		PlantedVariants: []PlantedVariant{planted},
		Reference:       refInput,
		Sequencer: SequencerManifest{
			Platform:      scenario.Sequencer.Platform,
			Coverage:      coverage,
			ReadLength:    readLength,
			PairedEnd:     paired,
			Bases:         bases,
			RegionFocused: fmt.Sprintf("%s:%d-%d", region.Chromosome, region.Start, region.End),
		},
		ComputeMetrics: ComputeMetrics{
			ReadsGenerated: readsGenerated,
			WallClockNS:    time.Since(start).Nanoseconds(),
		},
	}
	if scenario.Trio != nil && scenario.Trio.Generate {
		manifest.Trio = &TrioManifest{
			MotherCarries: false,
			FatherCarries: false,
			DeNovo:        scenario.Trio.DeNovo,
		}
		if simResult != nil {
			manifest.Trio.Population = string(simResult.Population)
			manifest.Trio.DemographicModel = simResult.Preset
			manifest.Trio.BackgroundVariants = simResult.NumVariants
		}
	}
	return manifest, nil
}

// scenarioRegion returns the 0-based half-open ±25 kb window around the
// variant. With a reference it is clipped to the contig, and the FASTA's
// name of the contig is returned too.
func (g *ScenarioGenerator) scenarioRegion(spec *PlantedVariantSpec) (*genomics.GenomicRegion, string, error) {
	varPos0 := spec.Pos - 1 // convert to 0-based
	region := &genomics.GenomicRegion{
		Chromosome: spec.Chrom,
		Start:      max(varPos0-windowHalfSize, 0),
		End:        varPos0 + windowHalfSize,
	}
	if g.referencePath == "" {
		return region, spec.Chrom, nil
	}
	fp, err := reference.NewFastaParser(g.referencePath)
	if err != nil {
		return nil, "", fmt.Errorf("reference %s: %w", g.referencePath, err)
	}
	index, err := fp.LoadIndex()
	if err != nil {
		return nil, "", fmt.Errorf("reference %s: %w", g.referencePath, err)
	}
	name, ok := matchChromosome(index, spec.Chrom)
	if !ok { // CheckPlantedRef has already rejected this
		return nil, "", fmt.Errorf("chromosome %q is not in reference %s", spec.Chrom, g.referencePath)
	}
	region.End = min(region.End, index.Chromosomes[name].Length)
	return region, name, nil
}

// writeReads runs the engine with cfg and opts and writes member's reads to
// out: paired reads go to out.r1 and out.r2, or interleaved to out.r1 when
// r2 is empty. It returns the number of FASTQ records written and records
// the files in rec. When regenerating (g.observe set) nothing is written:
// the bases of every read go to g.observe.
func (g *ScenarioGenerator) writeReads(cfg types.SequencerConfig, opts genomics.GenerateOpts, member string, out output, rec *runRecord) (int, error) {
	if out.r2 != "" && !cfg.PairedEnd {
		return 0, fmt.Errorf("an R2 file (%s) needs paired-end reads", out.r2)
	}
	if out.r2 != "" && strings.EqualFold(cfg.Output.Compression, "gzip") {
		return 0, fmt.Errorf("separate R1/R2 scenario files need uncompressed output (output compression %q)", cfg.Output.Compression)
	}
	engine := g.newEngine(cfg)
	defer engine.Close()
	if g.observe != nil {
		obs := &seqObserver{member: member, observe: g.observe}
		if err := engine.GenerateToWriter(obs, opts); err != nil {
			return 0, err
		}
		return obs.lines / 4, nil
	}

	paths := []string{out.r1}
	if out.r2 != "" {
		paths = append(paths, out.r2)
	}
	files := make([]*os.File, 0, len(paths))
	writers := make([]*bufio.Writer, 0, len(paths))
	defer func() {
		for _, f := range files {
			f.Close() // #nosec G104 -- closed explicitly on success
		}
	}()
	for _, p := range paths {
		f, err := os.Create(p) // #nosec G304 -- path from the scenario spec
		if err != nil {
			return 0, fmt.Errorf("creating %s: %w", p, err)
		}
		files = append(files, f)
		writers = append(writers, bufio.NewWriterSize(f, 256*1024))
	}
	split := &recordSplitter{out: make([]io.Writer, len(writers))}
	for i, w := range writers {
		split.out[i] = w
	}
	if err := engine.GenerateToWriter(split, opts); err != nil {
		return 0, err
	}
	for i, w := range writers {
		if err := w.Flush(); err != nil {
			return 0, err
		}
		if err := files[i].Close(); err != nil {
			return 0, err
		}
	}
	files = nil
	if strings.EqualFold(cfg.Output.Compression, "gzip") {
		n := cfg.ReadsPerBatch * len(engineMates(cfg))
		rec.addFile(member, out.r1, int64(n), 0)
		return n, nil
	}
	for i, p := range paths {
		rec.addFile(member, p, split.records[i], split.bases[i])
	}
	return split.lines / 4, nil
}

// engineMates returns the mates of one read unit.
func engineMates(cfg types.SequencerConfig) []int {
	if cfg.PairedEnd {
		return []int{1, 2}
	}
	return []int{1}
}

// recordSplitter sends successive 4-line FASTQ records to its writers in
// turn (R1, R2, R1, ... for interleaved pairs; everything to one writer when
// it has one), and counts the lines.
type recordSplitter struct {
	out   []io.Writer
	lines int
	// records and bases count each writer's records and bases.
	records, bases [2]int64
}

func (s *recordSplitter) Write(p []byte) (int, error) {
	n := 0
	for len(p) > 0 {
		chunk := p
		i := bytes.IndexByte(p, '\n')
		if i >= 0 {
			chunk = p[:i+1]
		}
		if _, err := s.out[(s.lines/4)%len(s.out)].Write(chunk); err != nil {
			return n, err
		}
		n += len(chunk)
		p = p[len(chunk):]
		k := (s.lines / 4) % len(s.out)
		if s.lines%4 == 1 {
			s.bases[k] += int64(len(chunk))
		}
		if i >= 0 {
			switch s.lines % 4 {
			case 0:
				s.records[k]++
			case 1:
				s.bases[k]-- // the newline
			}
			s.lines++
		}
	}
	return n, nil
}

// trioPlan is everything a trio member's reads are made from.
type trioPlan struct {
	sim  *population.HaplotypeSimulationResult
	trio *TrioHaplotypes
	// base holds the window's bases [start, start+len(base)) with every
	// background site's REF and the de novo REF in place.
	base  []byte
	start int64
	// runSeed is the scenario seed; each member's reads use a seed derived
	// from it (memberSeed).
	runSeed uint64
}

// planTrio simulates the six haplotypes for the window, the trio's
// transmission, and the window's bases. With a reference, background
// variants are anchored to it (see anchorToReference); without one, the
// window is one synthetic sequence drawn from the engine (uniform bases, or
// the HMM model when configured) with the simulated REF alleles written in.
func (g *ScenarioGenerator) planTrio(cfg types.SequencerConfig, region *genomics.GenomicRegion, refContig string, dn *genomics.VariantRecord) (*trioPlan, error) {
	regionLen := region.End - region.Start
	// Six phased haplotypes (mother, father, two background) from the
	// Go structured coalescent, one locus spanning the read region.
	sim, err := population.SimulateHaplotypes(
		population.Population(cfg.Population),
		region.Chromosome,
		int(region.Start),
		int(regionLen),
		6,
		cfg.Seed^trioSimSeedSalt,
	)
	if err != nil {
		return nil, fmt.Errorf("trio scenario: %w", err)
	}
	// #nosec G404 G115 -- deterministic simulation from the run seed, not security
	rng := mathrand.New(mathrand.NewSource(int64(cfg.Seed)))
	trio, err := GenerateTrioHaplotypes(sim, dn, rng)
	if err != nil {
		return nil, fmt.Errorf("generating trio haplotypes: %w", err)
	}
	tp := &trioPlan{sim: sim, trio: trio, start: region.Start, runSeed: cfg.Seed}

	if g.referencePath != "" {
		rc := cfg
		rc.ReferencePath = g.referencePath
		rc.ReferenceRegion = fmt.Sprintf("%s:%d-%d", refContig, region.Start+1, region.End)
		e := genomics.NewEngineLog(rc, io.Discard)
		tp.base, err = e.ReferenceBases(refContig, region.Start, region.End)
		e.Close()
		if err != nil {
			return nil, fmt.Errorf("trio scenario: %w", err)
		}
		anchorToReference(trio, tp.base, region.Start)
		return tp, nil
	}

	e := g.newEngine(cfg)
	defer e.Close()
	seq, _, err := e.GenerateRawSequenceWith(e.Rand(backgroundStream), int(regionLen))
	if err != nil {
		return nil, fmt.Errorf("trio scenario: generating the window: %w", err)
	}
	if int64(len(seq)) != regionLen {
		return nil, fmt.Errorf("trio scenario: generated %d window bases, want %d", len(seq), regionLen)
	}
	tp.base = bytes.ToUpper(seq)
	writeSyntheticREF(trio, tp.base, region.Start)
	return tp, nil
}

// anchorToReference places the simulated background variants on the
// reference bases base (the window starting at 0-based start). The
// coalescent gives each site's position and which haplotypes carry the
// derived allele; with a reference, the ancestral (REF) allele is the
// reference base there, and the derived (ALT) allele is the simulated
// derived base, or the simulated ancestral base when the simulated derived
// base equals the reference (the two simulated bases differ, so ALT always
// differs from REF). Sites on a base other than A/C/G/T, multi-base sites
// whose REF does not match the reference, and sites outside the window are
// not used (their alleles are cleared).
func anchorToReference(trio *TrioHaplotypes, base []byte, start int64) {
	for i := range trio.Variants {
		v := &trio.Variants[i]
		off := int64(v.Position) - start
		if len(v.Alleles) < 2 || len(v.Alleles[0]) == 0 || len(v.Alleles[1]) == 0 ||
			off < 0 || off+int64(len(v.Alleles[0])) > int64(len(base)) {
			v.Alleles = nil
			continue
		}
		if len(v.Alleles[0]) != 1 || len(v.Alleles[1]) != 1 {
			if !strings.EqualFold(string(base[off:off+int64(len(v.Alleles[0]))]), v.Alleles[0]) {
				v.Alleles = nil
			}
			continue
		}
		b := base[off]
		if !isACGT(b) {
			v.Alleles = nil
			continue
		}
		anc, der := strings.ToUpper(v.Alleles[0]), strings.ToUpper(v.Alleles[1])
		if der[0] == b {
			der = anc
		}
		v.Alleles = []string{string(b), der}
		v.AncestralState = string(b)
	}
}

// writeSyntheticREF writes each background site's simulated REF allele, and
// the de novo REF, into the synthetic window base at their positions.
func writeSyntheticREF(trio *TrioHaplotypes, base []byte, start int64) {
	put := func(pos0 int64, ref string) {
		off := pos0 - start
		if off >= 0 && off+int64(len(ref)) <= int64(len(base)) {
			copy(base[off:], strings.ToUpper(ref))
		}
	}
	for _, v := range trio.Variants {
		if len(v.Alleles) >= 2 && len(v.Alleles[0]) > 0 && len(v.Alleles[1]) > 0 {
			put(int64(v.Position), v.Alleles[0])
		}
	}
	if dn := trio.DeNovoVariant; dn != nil && dn.Position > 0 && dn.Ref != "" && dn.Alt != "" {
		put(dn.Position-1, dn.Ref)
	}
}

func isACGT(b byte) bool { return b == 'A' || b == 'C' || b == 'G' || b == 'T' }

// memberHaplotypes returns member's two haplotype sequences: the window
// bases with the member's allele written at every background site, and,
// on the proband's first haplotype, the de novo ALT. Sites are applied
// right to left so indels do not shift the sites to their left; a site
// overlapping one already applied is skipped.
func (tp *trioPlan) memberHaplotypes(member string) ([2][]byte, error) {
	t := tp.trio
	var haps [2][]int
	switch member {
	case "proband":
		haps = [2][]int{t.ProbandH1, t.ProbandH2}
	case "mother":
		haps = [2][]int{t.MotherH1, t.MotherH2}
	case "father":
		haps = [2][]int{t.FatherH1, t.FatherH2}
	default:
		return [2][]byte{}, fmt.Errorf("unknown member %q", member)
	}

	type site struct {
		off      int64
		ref, alt []byte
		hapIdx   int // index into the haplotype; -1 for the de novo variant
	}
	var sites []site
	for vi, v := range t.Variants {
		if len(v.Alleles) < 2 || len(v.Alleles[0]) == 0 || len(v.Alleles[1]) == 0 {
			continue
		}
		sites = append(sites, site{int64(v.Position) - tp.start, []byte(v.Alleles[0]), []byte(v.Alleles[1]), vi})
	}
	if dn := t.DeNovoVariant; dn != nil && dn.Position > 0 && dn.Ref != "" && dn.Alt != "" {
		sites = append(sites, site{dn.Position - 1 - tp.start, []byte(dn.Ref), []byte(dn.Alt), -1})
	}
	sort.SliceStable(sites, func(i, j int) bool { return sites[i].off > sites[j].off })

	var out [2][]byte
	for h := range out {
		seq := append([]byte(nil), tp.base...)
		limit := int64(len(seq))
		for _, st := range sites {
			end := st.off + int64(len(st.ref))
			if st.off < 0 || end > limit {
				continue
			}
			if !bytes.EqualFold(seq[st.off:end], st.ref) {
				return [2][]byte{}, fmt.Errorf("%s haplotype %d: REF %s at window offset %d does not match the window base %s",
					member, h+1, st.ref, st.off, seq[st.off:end])
			}
			var alt bool
			if st.hapIdx >= 0 {
				alt = st.hapIdx < len(haps[h]) && haps[h][st.hapIdx] == 1
			} else {
				alt = member == "proband" && h == 0
			}
			if alt {
				seq, _ = seqread.Splice(seq, nil, int(st.off), len(st.ref), st.alt)
			}
			limit = st.off
		}
		out[h] = seq
	}
	return out, nil
}

// writeMemberReads samples member's reads from its two haplotypes with the
// engine (reads per haplotype in proportion to its length, so about half
// each) and writes them to out. It returns the number of FASTQ records.
func (g *ScenarioGenerator) writeMemberReads(cfg types.SequencerConfig, tp *trioPlan, member string, out output, rec *runRecord) (int, error) {
	haps, err := tp.memberHaplotypes(member)
	if err != nil {
		return 0, err
	}
	dir, err := os.MkdirTemp("", "phantom-scenario-")
	if err != nil {
		return 0, err
	}
	defer os.RemoveAll(dir)
	fa := filepath.Join(dir, member+".haplotypes.fa")
	if err := writeFASTA(fa, []string{member + "_h1", member + "_h2"}, haps[:]); err != nil {
		return 0, err
	}
	mc := cfg
	mc.ReferencePath, mc.ReferenceRegion, mc.StrictRef = fa, "", false
	// Bases come from the haplotypes; a model has already shaped the
	// synthetic window.
	mc.HMMModelPath, mc.ModelDir = "", ""
	mc.Seed = memberSeed(tp.runSeed, member)
	return g.writeReads(mc, genomics.GenerateOpts{}, member, out, rec)
}

// writeFASTA writes the sequences as a FASTA file with 60-base lines.
func writeFASTA(path string, names []string, seqs [][]byte) error {
	f, err := os.Create(path) // #nosec G304 -- path in our own temporary directory
	if err != nil {
		return err
	}
	w := bufio.NewWriter(f)
	for i, s := range seqs {
		fmt.Fprintf(w, ">%s\n", names[i])
		for a := 0; a < len(s); a += 60 {
			w.Write(s[a:min(a+60, len(s))]) // #nosec G104 -- bufio errors surface at Flush
			w.WriteByte('\n')               // #nosec G104
		}
	}
	if err := w.Flush(); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// memberSeed returns the seed of a trio member's reads: the run seed for the
// proband, and a SplitMix64 mix of it for the parents, so members draw
// independent fragment positions, qualities and errors.
func memberSeed(seed uint64, member string) uint64 {
	var salt uint64
	switch member {
	case "mother":
		salt = 0x4d4f544845520000 // "MOTHER"
	case "father":
		salt = 0x4641544845520000 // "FATHER"
	default:
		return seed
	}
	z := seed ^ salt
	z += 0x9e3779b97f4a7c15
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	z ^= z >> 31
	if z == 0 {
		z = 1 // 0 means "draw a seed"
	}
	return z
}
