//go:build !release

package clinical

import (
	"bufio"
	"context"
	"fmt"
	"io"
	mathrand "math/rand"
	"os"
	"sort"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
	"github.com/scttfrdmn/phantom-science/internal/population"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

const windowHalfSize = 25000 // 50 kb window around variant (±25 kb)

// trioSimSeedSalt decorrelates the haplotype simulation's RNG from the
// engines' read streams, which use the same run seed.
const trioSimSeedSalt = 0x5452494f48415053 // "TRIOHAPS"

// ScenarioGenerator orchestrates haplotype-aware clinical scenario generation.
type ScenarioGenerator struct {
	config        *types.SequencerConfig
	referencePath string
}

// NewScenarioGenerator creates a ScenarioGenerator.
func NewScenarioGenerator(config *types.SequencerConfig) *ScenarioGenerator {
	return &ScenarioGenerator{config: config}
}

// WithReference makes Run check the planted variant's REF allele against
// the reference FASTA at path (see CheckPlantedRef) before writing any
// reads; a mismatch fails the run. Read bases remain synthetic.
func (g *ScenarioGenerator) WithReference(path string) *ScenarioGenerator {
	g.referencePath = path
	return g
}

// Run executes a clinical scenario and returns the ground truth manifest.
func (g *ScenarioGenerator) Run(ctx context.Context, scenario *ClinicalScenario) (*GroundTruthManifest, error) {
	start := time.Now()

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

	// 2. Build GenomicRegion: ±25 kb window around variant.
	varPos0 := spec.Pos - 1 // convert to 0-based
	region := &genomics.GenomicRegion{
		Chromosome: spec.Chrom,
		Start:      varPos0 - windowHalfSize,
		End:        varPos0 + windowHalfSize,
	}
	if region.Start < 0 {
		region.Start = 0
	}

	// Override from scenario sequencer spec if provided.
	readLength := g.config.ReadLength
	if scenario.Sequencer.ReadLength > 0 {
		readLength = scenario.Sequencer.ReadLength
	}
	coverage := scenario.Sequencer.Coverage
	if coverage <= 0 {
		coverage = 30
	}

	regionLen := region.End - region.Start
	totalReads := int(int64(coverage) * regionLen / int64(readLength))
	if totalReads < 1 {
		totalReads = 1
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

	// Clone config for generation.
	cfg := *g.config
	if readLength > 0 {
		cfg.ReadLength = readLength
	}
	cfg.ReadsPerBatch = totalReads
	cfg.TotalBatches = 1

	// Resolve the seed once so the haplotype simulation and every trio
	// member's engine derive from it; the scenario is reproducible with it.
	runSeed, drawn, err := genomics.ResolveSeed(cfg.Seed)
	if err != nil {
		return nil, err
	}
	if drawn {
		fmt.Fprintf(os.Stderr, "   [info] seed %d (set seed: %d in the sequencer config to reproduce this scenario)\n", runSeed, runSeed)
	}
	cfg.Seed = runSeed

	var readsGenerated int
	var simResult *population.HaplotypeSimulationResult

	if scenario.Trio != nil && scenario.Trio.Generate && scenario.Trio.DeNovo {
		// 4. Trio path: haplotype-aware generation.
		// Six phased haplotypes (mother, father, two background) from the
		// Go structured coalescent, one locus spanning the read region.
		result, err := population.SimulateHaplotypes(
			population.Population(cfg.Population),
			spec.Chrom,
			int(region.Start),
			int(regionLen),
			6,
			runSeed^trioSimSeedSalt,
		)
		if err != nil {
			return nil, fmt.Errorf("trio scenario: %w", err)
		}
		simResult = result

		// #nosec G404 G115 -- deterministic simulation from the run seed, not security
		rng := mathrand.New(mathrand.NewSource(int64(runSeed)))
		trio, err := GenerateTrioHaplotypes(result, &varRecord, rng)
		if err != nil {
			return nil, fmt.Errorf("generating trio haplotypes: %w", err)
		}

		// Proband FASTQ.
		pReads, err := g.writeHaplotypeFASTQ(ctx, scenario.Output.ProbandFASTQ, &cfg, trio, "proband", region, totalReads)
		if err != nil {
			return nil, fmt.Errorf("writing proband FASTQ: %w", err)
		}
		readsGenerated += pReads

		// Mother FASTQ (no de novo).
		if scenario.Output.MotherFASTQ != "" {
			mReads, err := g.writeHaplotypeFASTQ(ctx, scenario.Output.MotherFASTQ, &cfg, trio, "mother", region, totalReads)
			if err != nil {
				return nil, fmt.Errorf("writing mother FASTQ: %w", err)
			}
			readsGenerated += mReads
		}

		// Father FASTQ (no de novo).
		if scenario.Output.FatherFASTQ != "" {
			fReads, err := g.writeHaplotypeFASTQ(ctx, scenario.Output.FatherFASTQ, &cfg, trio, "father", region, totalReads)
			if err != nil {
				return nil, fmt.Errorf("writing father FASTQ: %w", err)
			}
			readsGenerated += fReads
		}
	} else {
		// 5. Simple proband FASTQ via Engine.GenerateToWriter + RefRegion.
		engine := genomics.NewEngine(cfg)
		defer engine.Close()
		opts := genomics.GenerateOpts{
			Variants:  []genomics.VariantRecord{varRecord},
			RefRegion: region,
		}
		f, err := os.Create(scenario.Output.ProbandFASTQ) // #nosec G304
		if err != nil {
			return nil, fmt.Errorf("creating proband FASTQ: %w", err)
		}
		bw := bufio.NewWriterSize(f, 256*1024)
		if err := engine.GenerateToWriter(bw, opts); err != nil {
			f.Close()
			return nil, fmt.Errorf("generating proband reads: %w", err)
		}
		if err := bw.Flush(); err != nil {
			f.Close()
			return nil, err
		}
		if err := f.Close(); err != nil {
			return nil, err
		}
		readsGenerated = totalReads
	}

	// 7. Build manifest.
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
	manifest := &GroundTruthManifest{
		GeneratedAt:     time.Now().UTC(),
		ProbandID:       g.config.SampleName,
		PlantedVariants: []PlantedVariant{planted},
		Reference:       g.referencePath,
		Sequencer: SequencerManifest{
			Platform:      scenario.Sequencer.Platform,
			Coverage:      coverage,
			ReadLength:    readLength,
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

// writeHaplotypeFASTQ generates a haplotype-aware FASTQ file for one trio member.
func (g *ScenarioGenerator) writeHaplotypeFASTQ(
	_ context.Context,
	path string,
	cfg *types.SequencerConfig,
	trio *TrioHaplotypes,
	member string,
	region *genomics.GenomicRegion,
	readCount int,
) (int, error) {
	engine := genomics.NewEngine(*cfg)
	defer engine.Close()

	f, err := os.Create(path) // #nosec G304
	if err != nil {
		return 0, err
	}
	defer f.Close()

	bw := bufio.NewWriterSize(f, 256*1024)
	written, err := generateHaplotypeFASTQ(bw, engine, cfg, trio, member, region, cfg.ReadLength, readCount)
	if err != nil {
		return 0, err
	}
	return written, bw.Flush()
}

// generateHaplotypeFASTQ writes haplotype-aware FASTQ reads to w.
// For each read it picks H1 or H2 (50/50) and injects alleles at their positions.
// The de novo variant is planted on H1 for the proband only.
func generateHaplotypeFASTQ(
	w io.Writer,
	engine *genomics.Engine,
	cfg *types.SequencerConfig,
	trio *TrioHaplotypes,
	member string,
	region *genomics.GenomicRegion,
	readLength int,
	readCount int,
) (int, error) {
	var h1, h2 []int
	switch member {
	case "proband":
		h1, h2 = trio.ProbandH1, trio.ProbandH2
	case "mother":
		h1, h2 = trio.MotherH1, trio.MotherH2
	case "father":
		h1, h2 = trio.FatherH1, trio.FatherH2
	default:
		return 0, fmt.Errorf("unknown member %q", member)
	}

	regionLen := region.End - region.Start - int64(readLength)
	if regionLen <= 0 {
		return 0, fmt.Errorf("region too small for read length %d", readLength)
	}

	// Every site that can be written into a read: background variants
	// (0-based positions) and the de novo variant (1-based → 0-based).
	type site struct {
		pos      int64
		ref, alt []byte
		hapIdx   int // index into the haplotype; -1 for the de novo variant
	}
	var sites []site
	for vi, v := range trio.Variants {
		if len(v.Alleles) < 2 || len(v.Alleles[0]) == 0 || len(v.Alleles[1]) == 0 {
			continue
		}
		sites = append(sites, site{int64(v.Position), []byte(v.Alleles[0]), []byte(v.Alleles[1]), vi})
	}
	if dn := trio.DeNovoVariant; dn != nil && dn.Chromosome != "" && dn.Position > 0 && len(dn.Ref) > 0 && len(dn.Alt) > 0 {
		sites = append(sites, site{dn.Position - 1, []byte(dn.Ref), []byte(dn.Alt), -1})
	}
	// Apply right to left so that indels do not shift the offsets of the
	// sites to their left; a site overlapping one already applied is skipped.
	sort.SliceStable(sites, func(i, j int) bool { return sites[i].pos > sites[j].pos })
	spans := make([]seqread.Span, len(sites))
	for i, st := range sites {
		spans[i] = seqread.Span{Pos: st.pos, RefLen: len(st.ref), AltLen: len(st.alt)}
	}

	rng := engine.Rand(memberStream(member))
	written := 0

	for i := 0; i < readCount; i++ {
		// Random read start within region.
		readStart := region.Start + rng.Int64N(regionLen)

		// Generate slack bases so that deletions are absorbed and the read
		// is trimmed back to readLength, never padded with 'N'.
		n := readLength + seqread.WindowSlack(readStart, readLength, spans)
		seq, quals, err := engine.GenerateRawSequenceWith(rng, n)
		if err != nil {
			return written, fmt.Errorf("generating sequence: %w", err)
		}
		if len(quals) != len(seq) {
			quals = make([]byte, len(seq))
			for qi := range quals {
				quals[qi] = 30
			}
		}

		// Pick H1 or H2 (50/50).
		useH1 := rng.IntN(2) == 0
		hap := h2
		if useH1 {
			hap = h1
		}

		// Write REF or ALT at every site in the read (so variant sites have
		// deterministic, haplotype-accurate bases). The de novo site is ALT
		// on proband H1 only: VAF 0.5 in the proband, 0 in the parents.
		window := int64(n)
		limit := window
		for _, st := range sites {
			if st.hapIdx >= len(hap) {
				continue
			}
			off := st.pos - readStart
			end := off + int64(len(st.ref))
			if off < 0 || end > window || end > limit {
				continue
			}
			allele := st.ref
			switch {
			case st.hapIdx >= 0 && hap[st.hapIdx] == 1:
				allele = st.alt
			case st.hapIdx < 0 && member == "proband" && useH1:
				allele = st.alt
			}
			seq, quals = seqread.Splice(seq, quals, int(off), len(st.ref), allele)
			limit = off
		}
		if len(seq) < readLength { // cannot happen: the slack covers every deletion
			return written, fmt.Errorf("read shorter than %d bases after variant injection", readLength)
		}
		seq, quals = seq[:readLength], quals[:readLength]

		// Encode quality scores (Phred+33, capped at 41).
		for qi, q := range quals {
			quals[qi] = '!' + min(q, 41)
		}

		readID := fmt.Sprintf("@%s:%d:HXXXXXX:1:1101:%d:%d 1:N:0:1",
			cfg.SampleName, cfg.Lane,
			int(readStart%100000), int((readStart/100000)%100000))
		if _, err := fmt.Fprintf(w, "%s\n%s\n+\n%s\n", readID, seq, quals); err != nil {
			return written, err
		}
		written++
	}
	return written, nil
}

// memberStream returns the RNG stream of a trio member's reads, so members
// generated by engines with the same seed draw independent reads.
func memberStream(member string) uint64 {
	switch member {
	case "mother":
		return 2
	case "father":
		return 3
	default:
		return 1
	}
}
