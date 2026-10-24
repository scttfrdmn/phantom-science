# Scientific Validation Methodology

This document describes the statistical and simulation models in Phantom
Science: what each computes, where it is implemented, which test checks it,
and what has not been validated. Every number quoted here comes from a test
or command named next to it.

Tests use the `!release` build tag. To run the ones cited below:

```bash
CGO_ENABLED=0 go test -tags '!release' -v ./internal/genomics/population/ \
    ./internal/genomics/authenticity/statistics/ ./internal/testing/ \
    ./internal/genomics/nanopore/errmodel/
```

---

## 1. Real versus synthetic: the goal and its status

The aim is synthetic data that a classifier cannot tell from real data. Two
distributions P (real) and Q (synthetic) are indistinguishable when

    TV(P, Q) = ½ · Σ |P(x) − Q(x)| = 0

where TV is the total variation distance. The best classifier then has an
accuracy of exactly 50%.

`cmd/run_e2e_tests` measures this with classifiers on three feature sets:

- **sequence:** 4-mer counts, Multinomial Naive Bayes (Section 4.5);
- **quality:** a per-read quality profile (quality histogram, mean quality
  of 10 cycle blocks, low-quality runs, quality changes between cycles, N
  fraction, mean quality per called base), L2 logistic regression;
- **combined:** a stacked classifier on the cross-fitted sequence log-odds
  and the quality features.

Each comparison trains on half of each read set and reports accuracy on the
other half. Every feature set is calibrated in the same run:

- **Null:** two independent draws from the engine; must lie within 0.5 ± 4σ.
- **GC control:** one draw GC-skewed to 0.70 (`testing.GCSkewed`); the
  sequence-aware sets must exceed 0.9.
- **Quality control:** the same seed with the other quality model (binned or
  linear); the quality-aware sets must exceed 0.9.
- **Real versus real:** with `-real-fastq FILE`, alternate halves of the real
  file; must lie within the null band.

With `-reference FASTA [-region]` the generated reads are sampled from that
genome, as `generate --reference` does.

**Measured against real reads** (2026-10-20, `run_e2e_tests -platforms
novaseq -reads 200000 -real-fastq` on the first 200,000 R1 reads of NovaSeq
6000 run DRR276626; 200,000 reads held out, band ±0.0045; real versus real
100,000 held out, ±0.0063):

| Features | Null | Real vs real | Uniform bases | `-reference` hg38, 0.11.0 | `-reference` hg38, base-dependent qualities (#35) |
|---|---|---|---|---|---|
| sequence (4-mer NB) | 0.500 | 0.498 | 0.946 | 0.518 | 0.519 |
| quality (LR) | 0.502 | 0.498 | 0.580 | 0.586 | **0.524** |
| combined | 0.502 | 0.499 | 0.951 | 0.586 | **0.527** |

Read 2 and pairs (`run_e2e_tests -real-fastq2`, paired `-reference`
output against real `DRR276626_2`): read-2 quality 0.543, combined 0.538;
the `pair` feature set (both mates' quality features and how they go
together) 0.593, with a pair-set null of 0.500.

**These numbers are realism for the fitted run only.** The model was fitted
to DRR276626 and measured against it. Measured against three held-out
NovaSeq 6000 runs from other centres and library preps
(`scripts/acceptance/illumina_heldout.sh`; ERR15600747, SureSelect QXT;
DRR394916, TruSeq PCR-free; DRR260182, TruSeq Nano; 200,000 pairs each,
measured 2026-10-23), compared with how far the real runs are from each
other (all six pairs of the four runs):

| Feature set | Generated vs held-out runs (min / median / max) | Real run vs real run (min / median / max) |
|---|---|---|
| read-1 sequence | 0.512 / 0.514 / 0.532 | 0.509 / 0.520 / 0.531 |
| read-1 quality | 0.619 / 0.672 / 0.762 | 0.573 / 0.646 / 0.767 |
| read-2 quality | 0.664 / 0.722 / 0.727 | 0.553 / 0.630 / 0.729 |
| pair | 0.718 / 0.770 / 0.807 | 0.618 / 0.677 / 0.780 |

Real-versus-real within each run (alternate halves) was 0.497–0.503 in all
28 comparisons.

- **Sequence** meets the realistic target: generated reads are as close to
  each run as real runs are to each other.
- **Qualities** behave like a copy of DRR276626. Generated reads are as far
  from each held-out run as DRR276626 itself is (read-1 quality 0.762 vs
  0.767, 0.672 vs 0.677, 0.619 vs 0.627). Real runs differ a great deal
  (0.57–0.77), and DRR276626 is the best of the four: mean Q 36.2 against
  34.6–35.6.
- **Pairs** were worse than DRR276626 on every run (+0.03 to +0.04): real
  mates' low-quality counts are correlated (0.28–0.53) and the model's mates
  were independent. The mates now share a quality class (a 6 × 6 joint class
  table fitted to DRR276626): count correlation −0.003 → 0.293 (real 0.282),
  both-mates-poor ratio 1.01 → 3.87 (real 3.95), and the pair classifier on
  DRR276626 0.593 → 0.551, with similar drops on every held-out run.

**Run profiles.** Since 0.13.0 the default profile is DRR394916, the most central against unfitted runs; the DRR276626 figures above apply to `--illumina-profile novaseq6000-drr276626`. Three runs are fitted with the same pipeline and shipped as
profiles; DRR260182, SRR30646135 (NIST GIAB) and ERR12737258 (MRC WIMM) are
never fitted and serve as held-out runs. Pair classifier, profile × run
(`results/heldout/profile_matrix.tsv`; real run vs real run 0.581–0.780,
median 0.675):

| Profile | DRR276626 | DRR394916 | ERR15600747 | DRR260182 | SRR30646135 | ERR12737258 |
|---|---|---|---|---|---|---|
| drr276626 | **0.551** | 0.742 | 0.777 | 0.678 | 0.752 | 0.690 |
| drr394916 | 0.757 | **0.570** | 0.654 | 0.655 | 0.640 | 0.648 |
| err15600747 | 0.796 | 0.662 | **0.583** | 0.698 | 0.703 | 0.707 |

Each profile is close to its own run, and against every other run it scores
about what its real run scores (within −0.003 to +0.037). Choosing the
profile nearest a held-out run in mean quality picks DRR394916 for all three,
scoring 0.640–0.655: within the range of the other real runs against those
runs, above their median. Limitation: read-2 G calls are poorer in DRR394916
and ERR15600747 than their profiles make them (mean Q 0.20 and 0.33 lower).

Results: `scripts/acceptance/calibrate/illumina/results/heldout/`
(`classifier_spread.tsv`, `classifier_matrix.tsv`,
`calibration_stats.tsv`).

In the same runs the GC control scored 0.973–0.994 and the quality control
1.000.

- **Uniform bases** are separated by composition: CpG 4-mers are about ten
  times more frequent than in human DNA.
- **Reads sampled from hg38** are close to real by 4-mer content (0.518).
  What remains: satellite-type 4-mers (GGAA, ATTC, TGGA) and CpG 4-mers are
  more frequent in real reads, and real GC is 0.414 against 0.408.
- **Qualities now depend on the read's bases (#35)**, which brought the
  quality classifier from 0.586 to 0.524. Real NovaSeq qualities depend on
  the base's context more than on the base itself (R1 P(Q25) 0.047 after TT,
  0.019 after CA; quality drops after GG), and the effect differs between
  good and poor reads, so the model multiplies each read class's Q11/Q25
  probabilities by factors of the 4-base context (two bases before, the
  base, one after), fitted to DRR276626.
- **What still separates qualities:** the arrangement of low calls along
  the read (#36; `scripts/acceptance/calibrate/illumina/results/heldout/structure_compare.tsv`).
  In all six runs:
  - good reads (1–4 calls below Q37) rarely repeat a Q11: P(Q11 | Q11)
    0.005–0.009, against 0.022–0.070 in the profiles;
  - Q11 is followed by Q25 more often, especially in cycles 1–30: 0.12–0.15,
    against 0.06–0.10;
  - good reads have their low calls at the read's end, while poor reads
    plateau after cycle 60; the model gives every read the same shape;
  - low calls cluster 2–20 cycles apart, and long low segments start
    mid-read (median cycle about 105, 9% reaching the end).

  Per-class cycle curves with a fitted carry between bins reproduce all of
  this except long runs. They improve each profile against its own run
  (read 2 −0.007 to −0.014, pair −0.008 to −0.013) but not against the
  never-fitted runs taken together (medians: read 1 +0.010, read 2 −0.006,
  pair +0.004), so they are not used. Each profile becomes closer to a copy
  of its own run, and real runs differ that much from each other. Part of
  the remaining per-read signal is in per-base mean quality and read-level
  effects (GC-rich, homopolymer and poly-G reads are poorer in real data).
  Read 2 also shows poly-G 4-mers in real data.
- A single region biases composition (chr22: sequence 0.618), so the
  whole-genome reference is the meaningful comparison.

The data is not yet indistinguishable from real reads, and these numbers
say by how much and why.

---

## 2. HMM sequence model

`phantom-science train` builds a model and `generate --hmm-model-path` draws
bases from it. Without a model, bases are uniform random.

### States and transitions

The hidden states are all (k−1)-mers over {A, C, G, T}. For k = 4 (the
default) there are 4³ = 64 states. T[s][b] = P(next base = b | state s), so
each state has 4 transition probabilities. Each state also has a 42-entry
emission vector for Phred qualities 0 to 41.

### Training

Training is a single counting pass, not iterative EM.

- `train` reads each VCF record's REF allele and first ALT allele and
  concatenates them. It does not read a reference genome, so the sequence
  context around variants is not learned.
- The adaptive model (`--model-type adaptive`, the default) streams that
  sequence in 100 KB chunks. Each chunk is prefixed with the last k−1 bases
  of the previous one, so chunked training equals training on the whole
  sequence. It counts each (k−1)-mer context followed by its next base.
- Both models estimate P(b | s) = (c(s,b) + 1) / (c(s) + 4), a pseudocount of
  1 per base, so unseen transitions keep a small probability. The start
  state is drawn from the smoothed (k−1)-mer frequencies. `markov_test.go`
  checks every cell against a direct tally.
- The k-mer model (`--model-type kmer`, k from 2 to 8) holds the sequence in
  memory and trains with `TrainFromReferenceGenome`.
- A VCF carries no base qualities, so emissions are not learned from it.

A model trained on little data mostly reproduces the smoothed counts. Before
0.10.0, missing normalisation and an off-by-one context made a model trained
on three short alleles generate reads of a single repeated base; now the same
input gives balanced bases (A 14890, C 14775, G 15186, T 15149 of 60,000).

### Generation-time reweighting in the adaptive model

`AdaptiveHMMModel` (`internal/genomics/authenticity/adaptive_hmm_model.go`)
still reweights transitions during generation for sequence complexity and
against homopolymer runs of 5 or more. These are modelling choices; no test
compares their effect with real data. The GC-content and CpG targets of
earlier versions were removed.

---

## 3. Population genetics: structured coalescent

`internal/genomics/population` simulates unlinked loci under a
multi-population coalescent with exponential growth, bottlenecks and
population splits (`SimulateJointPanel`). Coalescence uses the diploid rate,
with E[T₂] = 2N generations for a pair of lineages
(`TestStructured_ConstantN_PairwiseCoalescence`).

### Presets

`PresetDemography` (`presets.go`) models four 1000 Genomes populations as one
tree, with no migration and no admixture:

| Population | Ne | Growth per generation (last 204 generations) | Bottleneck | Split |
|---|---|---|---|---|
| YRI | 12,300 | 1.66% | 1000 for 20 generations from 2000 | root |
| CEU | 11,418 | 0.38% | 2000 for 400 generations from 1600 | from YRI, 2800 generations ago |
| CHB | 9,025 | 0.48% | 1500 for 100 generations from 1800 | from CEU, 1800 generations ago |
| MXL | 8,500 | 0.45% | 1200 for 100 generations from 1500 | from CHB, 600 generations ago |

Mutation rate 1.25 × 10⁻⁸ per bp per generation. Loci are simulated without
recombination, so the coalescent produces no linkage disequilibrium within a
locus. Other codes map to the closest preset as proxies: LWK, GWD and AFR to
YRI; GBR, TSI and EUR to CEU; JPT, CHS, EAS and ASN to CHB; PUR and AMR to
MXL. MXL is a Native-American-like branch, and its European and African
admixture is not modelled. `presets.go` documents where each parameter comes
from and which were calibrated.

### What the tests check

`presets_test.go` simulates YRI, CEU and CHB together (50 diploids each, 2000
loci of 10 kb). From one run on 2026-10-18:

| Check | Test | Result |
|---|---|---|
| π ordering YRI > CEU, CHB (at 4σ) | `TestPresets_DiversityOrdering` | π = 6.04, 5.62, 5.81 × 10⁻⁴ per bp |
| WC F_ST YRI–CEU in 0.10–0.16 | `TestPresets_FSTLiteratureRanges` | 0.143 ± 0.004 |
| WC F_ST YRI–CHB in 0.09–0.15 | same | 0.136 ± 0.003 |
| WC F_ST CEU–CHB in 0.06–0.12 | same | 0.096 ± 0.002 |
| Two samples of one population: F_ST < 0.02 | `TestPresets_SamePopulationNullFST` | 0.000 |

The literature ranges are the ones recorded in
`internal/population/fst_analysis.py` (removed after v0.9.0). Split times and bottleneck durations
were calibrated to fall inside them, so these checks confirm the calibration
rather than test it independently.

### Where it is used

`NewPopulationPanel` and `NewJointPanels` build variant panels (sites,
frequencies, genotypes) from one seeded simulation. `generate` and `train` do
not use the coalescent. Their `--population` flag only names a model file.
Trio generation in `scenario` uses it: `internal/population` runs
`SimulateJointPanel` with the scenario population's preset over one locus
spanning the scenario region and reads each haplotype's derived alleles from
`JointSite.Carriers`.

---

## 4. Statistics

### 4.1 Jensen-Shannon divergence

    JS(P ‖ Q) = ½ KL(P ‖ M) + ½ KL(Q ‖ M)    where M = (P + Q) / 2
    KL(P ‖ Q) = Σ P(x) log₂(P(x) / Q(x))

Range [0, 1] in bits (the default base is 2).

Implementation: `internal/genomics/authenticity/statistics/jensen_shannon.go`.
Tests: `jensen_shannon_test.go`.

### 4.2 F_ST (Weir and Cockerham 1984)

Inputs: allele frequency pⱼ and diploid sample size nⱼ for each of r
populations.

    N   = Σ nⱼ                                  total individuals
    p̄   = Σ nⱼ·pⱼ / N                           weighted mean frequency
    nᶜ  = (N − Σ nⱼ²/N) / (r−1)                 size-correction factor

    MSG = Σ nⱼ·(pⱼ − p̄)² / (r−1)               mean square between populations
    MSW = Σ nⱼ·pⱼ·(1−pⱼ) / (N−r)               mean square within populations

    F_ST = (MSG − MSW) / (MSG + (nᶜ−1)·MSW)

The result is clamped to [0, 1]. `Valid` is false when MSG ≤ 0 (no
between-population variance). `ComputeMultiLocusFST` averages valid loci.

Implementation: `statistics/fst.go`. Tests: `fst_test.go`.

### 4.3 Hardy-Weinberg equilibrium test

Inputs: genotype counts n_AA, n_Aa, n_aa (N individuals).

    p̂    = (2·n_AA + n_Aa) / (2·N)
    E_AA = N·p̂²,  E_Aa = 2·N·p̂·(1−p̂),  E_aa = N·(1−p̂)²

    χ² = Σ (O − E)² / E    (1 degree of freedom)
    p  = erfc(√(χ²/2))     (exact for 1 df)

Implementation: `statistics/hwe.go`. Tests: `hwe_test.go`.

### 4.4 Linkage disequilibrium (r², D, D′)

Inputs: phased haplotype counts n_AB, n_Ab, n_aB, n_ab (N haplotypes).

    p_A = (n_AB + n_Ab) / N,    q_a = 1 − p_A
    p_B = (n_AB + n_aB) / N,    q_b = 1 − p_B

    D    = n_AB/N − p_A·p_B
    r²   = D² / (p_A·q_a·p_B·q_b)
    D_max = min(p_A·q_b, q_a·p_B)    if D > 0
          = max(−p_A·p_B, −q_a·q_b)  if D ≤ 0
    D′   = D / D_max                  (range [−1, 1])

Implementation: `statistics/ld.go`. Tests: `ld_test.go`.

### 4.5 Multinomial Naive Bayes classifier

Features: a 256-element 4-mer count vector. A sequence of length L gives L−3
counts.

    θ_{c,k} = (count_k(c) + α) / (Σⱼ count_j(c) + α·256)    with α = 1
    class   = argmax_c  [ log P(c) + Σₖ xₖ · log θ_{c,k} ]

`CrossValidate` runs k-fold stratified cross-validation and returns the mean
and standard deviation of the per-fold accuracy.

Implementation: `internal/testing/naive_bayes.go`. Tests:
`naive_bayes_test.go`.

### 4.6 1-Wasserstein distance (supplementary)

`CompareDistributions` reports it next to the JS divergence. It integrates
|CDF_P − CDF_Q| over the N k-mers sorted lexicographically and divides by
N − 1. This is a Wasserstein distance with uniform cost on lexicographic rank,
not an earth-mover distance with Hamming cost between k-mers. Use the JS
divergence for biological comparisons.

---

## 5. Nanopore models

### 5.1 Basecaller error and quality model

`internal/genomics/nanopore/errmodel` (since 0.8.0) applies per-base errors
in one pass over the true sequence. Each base is deleted with probability
Del + homopolymer excess, substituted with probability Sub, and followed by an
inserted base with probability Ins:

| Tier | Sub | Ins | Del | Homopolymer excess cap | Modal read Q |
|---|---|---|---|---|---|
| sup | 0.8% | 0.8% | 1.0% | 4% | 20 |
| hac | 1.6% | 1.4% | 2.0% | 6% | 15 |
| fast | 5.0% | 3.5% | 6.5% | 10% | 10 |

The homopolymer excess grows with run length L as cap · (1 − exp(−(L−1)/4)).
Read quality is the ONT mean Q, −10·log10(mean(10^(−q/10))). These values are
taken from the project roadmap and published ranges (sources are in the
package comment). They are not fitted to real data, and no test compares
the output with real basecalls. `errmodel_test.go` checks that the realised
rates, homopolymer deletions and qualities match the parameters
(`TestRatesPerTier`, `TestHomopolymerDeletions`, `TestQualities`).

### 5.2 Read lengths

`internal/genomics/readlen` samples fixed, log-normal (from N50 and median,
or mean and SD), gamma or empirical lengths, with optional truncation.
After generation, reads are checked against the configured distribution
(bounds and a 4σ test on the mean). A standalone `validate` reports length
statistics (count, min, max, mean, N50), since it does not know the intended
distribution.

### 5.3 Raw signal

`internal/genomics/signal` reimplements squigulator's algorithm: k-mer
pore-model levels, gamma dwell, level and white noise, baseline drift,
stalls, skips, an adapter leader, poly-A/T plateaus and int16 digitisation.
The pore models are ONT's published tables, scaled to pA with the fixed
standardisation in Dorado's `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` config
(mean 93.6924 pA, stdev 23.5067 pA).

`scripts/acceptance/dorado_identity.sh` basecalls the signal with Dorado
SUP and aligns the calls to the source sequences. Measured on 2026-10-18
with Dorado 2.1.2:

| Signal | Median identity |
|---|---|
| Real R10.4.1 reads (`ont-open-data/pgx_as_2025.07`, against GRCh38) | 98.52% |
| squigulator 0.5.0 | 91.47% |
| `signal`, default DNA profile | 90.79% |
| `signal`, the 0.8.0 DNA profile | 88.05% |

The DNA profile's dwell, stall, skip and noise values were chosen in these
runs: short dwells and skips cost accuracy, and the real dwell tail is
reproduced by stalls on a low-variance gamma dwell (mean dwell 13.0
samples/base, CV 0.89; real 13.1, CV 1.08). `dorado_identity.sh` therefore
passes signal within one point of squigulator on the same truth rather than
at a fixed identity.

**Where the gap comes from** (#28; `scripts/acceptance/calibrate/`, 2000
reads per condition, bootstrap 95% CIs; results in
`scripts/acceptance/calibrate/results/dna.tsv`, measured 2026-10-20):

| Condition | Median identity [95% CI] |
|---|---|
| Real reads (4120 aligned to GRCh38) | 98.54 [98.52, 98.58] |
| `signal`, uniform random input | 90.80 [90.76, 90.82] |
| `signal`, GRCh38 pseudo-reads | 91.75 [91.67, 91.82] |
| squigulator, GRCh38 | 93.74 [93.68, 93.80] |
| squigulator with ONT's `9mer_levels_v1`, GRCh38 | 91.63 [91.52, 91.76] |
| `signal` with squigulator's levels, GRCh38 | 93.81 [93.74, 93.87] |
| `level_tau` 0.5 / 1 / 2 / 4 | 91.49 / 91.29 / 87.77 / 72.90 |

- **Composition, +1.0 point:** real-genome input against uniform random
  bases.
- **The level table, +2.1 points:** squigulator's R10 levels instead of
  ONT's `9mer_levels_v1`; squigulator loses the same with ONT's table. This
  was all of squigulator's lead over `signal`.
- **Ruled out:** level smoothing between k-mers (`level_tau`, worse at every
  τ, so it is 0 by default), per-sample noise matched to real reads, a level
  shift, and a per-k-mer stdv.
- **Remaining ~4.6 points:** the k-mer level model itself, which squigulator
  shares. Real per-base levels fit either table at r 0.70–0.72, against 0.85
  for synthetic signal; CpG k-mers sit 3.8 pA lower in native (methylated)
  DNA. Closing this needs a level table fitted by a resquiggler; refitting
  ONT's table from Dorado's move boundaries made identity worse (77.5%).

**RNA004** (`rna004_130bps_4khz`, Dorado `rna004_130bps_sup@v5.3.0`, the six
100-nt control constructs of `ont-open-data/rna-modbase-validation_2025.03`,
`results/rna.tsv`): the profile is calibrated against the real reads
(adapter 0.7 s at 62 pA, scaling 17/80 pA from Dorado's model config, 3.0 pA
noise through `model_stdv`, no white noise, drift or skips). Real reads
reach 97.67% median identity; calibrated `signal` 72.92% [72.09, 73.75]
(72% of reads align), up from 63.41% (62%); squigulator 82.86%. The main
missing piece is a per-k-mer stdv, which ONT's RNA004 table does not carry:
with squigulator's per-k-mer-stdv table, 97% of reads align and identity is
77.8%. Synthetic RNA004 signal is not yet basecalled like real signal.

### 5.4 cDNA transcriptome

`scripts/acceptance/transcriptome_bambu.sh` runs Pychopper, minimap2 and
Bambu on simulated cDNA reads and scores them against the simulation's
truth. On GRCh38 chr22 with GENCODE v49 (200,000 reads, novel rate 0.2,
2026-10-18): Bambu transcript counts correlate with the true read counts at
Pearson 0.872 (log1p) and Spearman 0.832; Bambu (NDR 1) reports novel
isoforms with precision 0.90 against the injected exon skips, and recovers
55% of the injected exon-skip isoforms that hold at least 5 reads and 5% of
their gene's reads.

---

## 6. Summary

| Model or metric | Implemented | Validated against | Not validated |
|---|---|---|---|
| JS divergence, F_ST, HWE, LD, Wasserstein | yes | unit tests of the formulas | — |
| Naive Bayes classifier | yes | null and GC-skew calibration | real vs synthetic accuracy |
| HMM sequence model | yes | unit tests | resemblance to real sequence |
| Structured coalescent presets | yes | literature F_ST ranges, π ordering (calibrated) | LD, admixture, other populations |
| Illumina qualities and errors | fixed quality ramp, no errors | — | everything |
| ONT error and quality model | yes | parameter round trip | real basecalls |
| Raw signal | yes | BLOW5 format (slow5tools, blue-crab); Dorado identity against squigulator; gap to real reads attributed (composition, level table, level model) | DNA identity ~6.8 points below real; RNA004 ~25 points below real |
| cDNA transcriptome | yes | Bambu counts and novel-isoform discovery | real cDNA libraries |

---

## 7. References

1. Weir, B.S. & Cockerham, C.C. (1984). Estimating F-statistics for the
   analysis of population structure. *Evolution* 38(6), 1358–1370.
2. Lewontin, R.C. (1964). The interaction of selection and linkage.
   *Genetics* 49(1), 49–67.
3. Hill, W.G. & Robertson, A. (1968). Linkage disequilibrium in finite
   populations. *Theoretical and Applied Genetics* 38(6), 226–231.
4. Bird, A.P. (1980). DNA methylation and the frequency of CpG in animal DNA.
   *Nucleic Acids Research* 8(7), 1499–1504.
5. Lander, E.S. et al. (2001). Initial sequencing and analysis of the human
   genome. *Nature* 409, 860–921.
6. Wigginton, J.E., Cutler, D.J. & Abecasis, G.R. (2005). A note on exact
   tests of Hardy-Weinberg equilibrium. *Am. J. Hum. Genet.* 76, 887–893.
7. McCallum, A. & Nigam, K. (1998). A comparison of event models for Naive
   Bayes text classification. *AAAI Workshop on Learning for Text
   Categorization.*
8. Tennessen, J.A. et al. (2012). Evolution and functional impact of rare
   coding variation from deep sequencing of human exomes. *Science* 337,
   64–69.
9. Gravel, S. et al. (2011). Demographic history and rare allele sharing
   among human populations. *PNAS* 108(29), 11983–11988.
10. Delahaye, C. & Nicolas, J. (2021). Sequencing DNA with nanopores:
    Troubles and biases. *PLoS One* 16, e0257521.
11. squigulator, the nanopore signal simulator whose algorithm
    `internal/genomics/signal` reimplements:
    https://github.com/hasindu2008/squigulator
