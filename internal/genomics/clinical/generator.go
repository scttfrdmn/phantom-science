//go:build !release

package clinical

import (
    "bufio"
    "context"
    "crypto/rand"
    "fmt"
    "io"
    mathrand "math/rand"
    "os"
    "time"

    "github.com/scttfrdmn/phantom-science/internal/genomics"
    "github.com/scttfrdmn/phantom-science/internal/genomics/diseases"
    "github.com/scttfrdmn/phantom-science/internal/population"
    "github.com/scttfrdmn/phantom-science/pkg/types"
)

const windowHalfSize = 25000 // 50 kb window around variant (±25 kb)

// ScenarioGenerator orchestrates haplotype-aware clinical scenario generation.
type ScenarioGenerator struct {
    config      *types.SequencerConfig
    projectRoot string
}

// NewScenarioGenerator creates a ScenarioGenerator.
func NewScenarioGenerator(config *types.SequencerConfig, projectRoot string) *ScenarioGenerator {
    return &ScenarioGenerator{config: config, projectRoot: projectRoot}
}

// Run executes a clinical scenario and returns the ground truth manifest.
func (g *ScenarioGenerator) Run(ctx context.Context, scenario *ClinicalScenario) (*GroundTruthManifest, error) {
    start := time.Now()

    // 1. Resolve variant.
    var spec *PlantedVariantSpec
    var gene string
    if scenario.Phenotype.PlantedVariant != nil {
        spec = scenario.Phenotype.PlantedVariant
        gene = scenario.Phenotype.ResolvedGene
    } else {
        cv := diseases.NewClinVarIntegrator("")
        var err error
        spec, gene, err = ResolveHPO(ctx, scenario.Phenotype.HPOTerms, cv)
        if err != nil {
            return nil, fmt.Errorf("resolving HPO terms: %w", err)
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

    var readsGenerated int

    if scenario.Trio != nil && scenario.Trio.Generate && scenario.Trio.DeNovo {
        // 4. Trio path: haplotype-aware generation.
        seed := int(time.Now().UnixNano() & 0x7fffffff)
        sim := population.NewHaplotypeSimulator(g.projectRoot)
        result, err := sim.SimulateHaplotypes(
            population.Population(g.config.Population),
            spec.Chrom,
            3, // 6 haplotypes
            int(regionLen),
            &seed,
            true,
        )
        if err != nil {
            return nil, fmt.Errorf("haplotype simulation: %w", err)
        }

        // #nosec G404 -- seed is used for deterministic simulation, not security
        rng := mathrand.New(mathrand.NewSource(int64(seed)))
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
    manifest := &GroundTruthManifest{
        GeneratedAt: time.Now().UTC(),
        ProbandID:   g.config.SampleName,
        PlantedVariants: []PlantedVariant{{
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
        }},
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

    qualityChars := "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
    coinBuf := make([]byte, 1)
    rsBuf := make([]byte, 8)
    written := 0

    for i := 0; i < readCount; i++ {
        // Random read start within region.
        rand.Read(rsBuf) // #nosec G104
        rsVal := uint64(rsBuf[0])<<56 | uint64(rsBuf[1])<<48 | uint64(rsBuf[2])<<40 |
            uint64(rsBuf[3])<<32 | uint64(rsBuf[4])<<24 | uint64(rsBuf[5])<<16 |
            uint64(rsBuf[6])<<8 | uint64(rsBuf[7])
        rsVal &= 0x7fffffffffffffff
        readStart := region.Start + int64(rsVal%uint64(regionLen))

        seq, hmmQuals, err := engine.GenerateRawSequence(readLength)
        if err != nil {
            return written, fmt.Errorf("generating sequence: %w", err)
        }

        // Pick H1 or H2 (50/50).
        rand.Read(coinBuf) // #nosec G104
        useH1 := coinBuf[0]&1 == 0
        var hap []int
        if useH1 {
            hap = h1
        } else {
            hap = h2
        }

        // Set alleles at background variant positions (always set REF or ALT
        // so that variant sites have deterministic, haplotype-accurate bases).
        for vi, v := range trio.Variants {
            if vi >= len(hap) {
                break
            }
            if len(v.Alleles) < 2 {
                continue
            }
            vPos0 := int64(v.Position) // HaplotypeVariant.Position is 0-based
            off := vPos0 - readStart
            if off < 0 || int(off)+len(v.Alleles[0]) > len(seq) {
                continue
            }
            ref := []byte(v.Alleles[0])
            alt := []byte(v.Alleles[1])
            if hap[vi] == 1 {
                seq = spliceAllele(seq, int(off), ref, alt) // inject ALT
            } else {
                seq = spliceAllele(seq, int(off), ref, ref) // enforce REF
            }
        }

        // Set base at de novo position: ALT on proband H1, REF everywhere else.
        // This ensures VAF = 0.5 for proband and 0 for mother/father.
        if trio.DeNovoVariant != nil {
            dn := trio.DeNovoVariant
            if dn.Chromosome != "" && dn.Position > 0 && len(dn.Ref) > 0 && len(dn.Alt) > 0 {
                off := (dn.Position - 1) - readStart // dn.Position is 1-based
                if off >= 0 && int(off)+len(dn.Ref) <= len(seq) {
                    if member == "proband" && useH1 {
                        seq = spliceAllele(seq, int(off), []byte(dn.Ref), []byte(dn.Alt))
                    } else {
                        seq = spliceAllele(seq, int(off), []byte(dn.Ref), []byte(dn.Ref))
                    }
                }
            }
        }

        // Encode quality scores.
        quals := make([]byte, readLength)
        if hmmQuals != nil {
            for qi, q := range hmmQuals {
                idx := int(q)
                if idx >= len(qualityChars) {
                    idx = len(qualityChars) - 1
                }
                quals[qi] = qualityChars[idx]
            }
        } else {
            for qi := range quals {
                quals[qi] = qualityChars[30]
            }
        }

        readID := fmt.Sprintf("@%s:%d:HXXXXXX:1:1101:%d:%d 1:N:0:1",
            cfg.SampleName, cfg.Lane,
            int(readStart%100000), int((readStart/100000)%100000))
        fmt.Fprintf(w, "%s\n%s\n+\n%s\n", readID, string(seq), string(quals))
        written++
    }
    return written, nil
}

// spliceAllele replaces bases at offset with alt, preserving original length.
func spliceAllele(seq []byte, offset int, ref, alt []byte) []byte {
    after := make([]byte, 0, len(seq))
    after = append(after, seq[:offset]...)
    after = append(after, alt...)
    after = append(after, seq[offset+len(ref):]...)
    if len(after) > len(seq) {
        after = after[:len(seq)]
    }
    for len(after) < len(seq) {
        after = append(after, 'N')
    }
    return after
}
