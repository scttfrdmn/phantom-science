# Scientific Validation Methodology

This document describes the mathematical foundations of every statistical metric
implemented in Phantom Science. It distinguishes between metrics that are fully
implemented, those that require real VCF data to produce a meaningful number, and
those that are documented approximations.

---

## 1. What is Statistical Indistinguishability?

Two distributions P (real genomic sequences) and Q (synthetic sequences) are
**statistically indistinguishable** if no classifier trained on samples from P and Q
can predict the source label with accuracy significantly above 50%. Formally:

    TV(P, Q) = ½ · Σ |P(x) − Q(x)| ≈ 0

where TV is the total variation distance. When TV(P, Q) = 0, P and Q are identical, and
the optimal classifier achieves exactly 50% accuracy (random guessing).

We operationalise this via a **Multinomial Naive Bayes classifier on 4-mer frequency
vectors** with 5-fold stratified cross-validation. If the measured accuracy is not
significantly above 50%, the sequences are considered indistinguishable at the k-mer level.

---

## 2. HMM Model Architecture

### States
All (k−1)-mers over {A, T, G, C} form the hidden states.
For k=4 (the default), there are 4³ = 64 states.

### Transitions
T\[s\]\[b\] = P(next base = b | current state = s), giving 4 transition probabilities per state.

### Emissions
E\[s\]\[q\] = P(Phred quality = q | state = s), for q ∈ {0, …, 41} — 42 values per state.

### Biological constraints enforced during generation

| Constraint | Value | Source |
|------------|-------|--------|
| GC content | 42% ± 1.5% | Lander et al. 2001 (human genome average) |
| CpG depletion | 25% of neutral expectation | Bird 1980 |
| Homopolymer limit | max run length 8 | Illumina accuracy ceiling |

### Training

Parameters are estimated in a **single pass** from k-mer frequencies observed in the
VCF REF + first ALT alleles, streamed in 100 KB chunks, with Laplace smoothing (α=0.1).
This is **not** iterative EM — the model is a direct frequency estimator.

---

## 3. Population-Specific Modeling

Allele frequencies are sourced from the **1000 Genomes Project Phase 3**
(Auton et al. 2015, *Nature* 526, 68–74) for 11 superpopulations:

| Code | Population |
|------|-----------|
| YRI | Yoruba in Ibadan, Nigeria |
| LWK | Luhya in Webuye, Kenya |
| GWD | Gambian in Western Divisions, The Gambia |
| CEU | Utah residents (CEPH) with N/W European ancestry |
| GBR | British in England and Scotland |
| TSI | Toscani in Italy |
| IBS | Iberian populations in Spain |
| CHB | Han Chinese in Beijing |
| JPT | Japanese in Tokyo |
| CHS | Southern Han Chinese |
| MXL | Mexican Ancestry in Los Angeles |

Selected variants with known population differentiation include:

- **APOE4** (rs429358): risk variant for Alzheimer's; elevated in YRI (~22%) vs CEU (~14%)
- **ACTN3** R577X (rs1815739): athletic performance; variable across populations
- **LCT** (rs4988235): lactase persistence; high in CEU (~77%), absent in CHB
- **CYP2C19** *2 (rs4244285): drug metabolism; elevated in CHB/JPT vs CEU

LD structure uses coalescent simulation via **msprime** (Kelleher et al. 2016,
*PLOS Computational Biology*).

---

## 4. Statistical Validation Suite

### 4.1 Jensen-Shannon Divergence

**Formula:**

    JS(P ‖ Q) = ½ KL(P ‖ M) + ½ KL(Q ‖ M)    where M = (P + Q) / 2

    KL(P ‖ Q) = Σ P(x) log₂(P(x) / Q(x))

**Range:** [0, 1] (in bits when using log₂).

**Implementation:** `internal/genomics/authenticity/statistics/jensen_shannon.go`

**Status:** ✅ Fully implemented and tested.

---

### 4.2 F_ST — Population Differentiation

Implements the **Weir-Cockerham (1984)** F_ST estimator.

**Inputs:** Allele frequency pⱼ and diploid sample size nⱼ for each of r populations.

**Derivation:**

    N   = Σ nⱼ                                  total individuals
    p̄   = Σ nⱼ·pⱼ / N                           weighted mean frequency
    nᶜ  = (N − Σ nⱼ²/N) / (r−1)                 size-correction factor

    MSG = Σ nⱼ·(pⱼ − p̄)² / (r−1)               mean square between populations
    MSW = Σ nⱼ·pⱼ·(1−pⱼ) / (N−r)               mean square within populations

    F_ST = (MSG − MSW) / (MSG + (nᶜ−1)·MSW)

Result clamped to [0, 1]. Returns `Valid=false` when MSG ≤ 0 (no between-population variance).

**Implementation:** `internal/genomics/authenticity/statistics/fst.go`

**Status:** ✅ Fully implemented and tested.

**Reference:** Weir, B.S. & Cockerham, C.C. (1984). Estimating F-statistics for the
analysis of population structure. *Evolution* 38(6), 1358–1370.

---

### 4.3 Hardy-Weinberg Equilibrium Test

**Inputs:** Observed genotype counts n_AA, n_Aa, n_aa (N = total individuals).

**Derivation:**

    p̂  = (2·n_AA + n_Aa) / (2·N)    MLE allele frequency
    E_AA = N·p̂²
    E_Aa = 2·N·p̂·(1−p̂)
    E_aa = N·(1−p̂)²

    χ² = (n_AA−E_AA)²/E_AA + (n_Aa−E_Aa)²/E_Aa + (n_aa−E_aa)²/E_aa    (df = 1)

    p-value = erfc(√(χ²/2))    [exact for df=1 via math.Erfc from Go stdlib]

**Implementation:** `internal/genomics/authenticity/statistics/hwe.go`

**Status:** ✅ Fully implemented and tested.

**Reference:** Wigginton, J.E., Cutler, D.J. & Abecasis, G.R. (2005). A note on exact
tests of Hardy-Weinberg equilibrium. *Am. J. Hum. Genet.* 76, 887–893.

---

### 4.4 Linkage Disequilibrium (r², D, D')

**Inputs:** Phased haplotype counts n_AB, n_Ab, n_aB, n_ab (N = total haplotypes).

**Derivation:**

    p_A = (n_AB + n_Ab) / N,    q_a = 1 − p_A
    p_B = (n_AB + n_aB) / N,    q_b = 1 − p_B
    p_AB = n_AB / N

    D    = p_AB − p_A·p_B
    r²   = D² / (p_A·q_a·p_B·q_b)

    D_max = min(p_A·q_b, q_a·p_B)       if D > 0
          = max(−p_A·p_B, −q_a·q_b)     if D ≤ 0
    D'   = D / D_max                     (range [−1, 1])

**Implementation:** `internal/genomics/authenticity/statistics/ld.go`

**Status:** ✅ Fully implemented and tested.

**References:**
- Lewontin, R.C. (1964). The interaction of selection and linkage. *Genetics* 49(1), 49–67.
- Hill, W.G. & Robertson, A. (1968). Linkage disequilibrium in finite populations.
  *Theoretical and Applied Genetics* 38(6), 226–231.

---

### 4.5 ML Classifier Accuracy (Multinomial Naive Bayes)

**Purpose:** Directly measure the ability of a classifier to distinguish real from
synthetic sequences. Accuracy ≈ 50% indicates indistinguishability.

**Feature extraction:** 256-element 4-mer frequency count vector. For a sequence
of length L, the vector has L−3 total counts (sliding window of width 4).

**Smoothed log-probabilities (Laplace, α=1):**

    θ_{c,k} = (count_k(c) + α) / (Σⱼ count_j(c) + α·256)

**Prediction:**

    log P(c | x) = log P(c) + Σₖ xₖ · log θ_{c,k}
    class = argmax_c  log P(c | x)

**Evaluation:** 5-fold stratified cross-validation. Returns mean ± std accuracy per fold.

**Implementation:** `internal/testing/naive_bayes.go`

**Status:** ✅ Implemented. Measured indistinguishability requires real VCF input.

---

### 4.6 1-Wasserstein Distance (Supplementary Metric)

Implemented in `CompareDistributions` as a supplementary metric alongside JS divergence.

**Note on normalisation:** The implementation integrates |CDF_P − CDF_Q| over N
lexicographically sorted support points and divides by N−1 (the correct maximum for
N ordered points). This is the 1-Wasserstein distance with **uniform transport cost
on lexicographic rank**. It is **not** the earth-mover distance with Hamming cost
between k-mer strings. For biological comparisons, JS divergence is preferred.

**Status:** ✅ Implemented (normalisation corrected in v0.7.0).

---

## 5. Current Validation Status

| Metric | Implementation | Requires Real Data | Notes |
|--------|---------------|-------------------|-------|
| k-mer JS divergence | ✅ | No (runs on generated data) | Mathematically correct |
| F_ST (Weir-Cockerham) | ✅ | Yes (real population allele frequencies) | Correct estimator |
| Hardy-Weinberg χ² | ✅ | Yes (real genotype counts) | Correct test |
| LD r²/D' | ✅ | Yes (real phased haplotypes) | Correct formula |
| ML accuracy (MNB, 5-fold CV) | ✅ | Yes (real sequences to test against) | ~50% = indistinguishable |
| Platform error models | Approximation | No | Reasonable; not empirically fitted |
| Population LD patterns | Approximation | No | Coalescent simulation; not matched to specific VCF |

---

## 6. How to Run Validation on Real Data

### Train a model from a VCF file

```bash
phantom-science train \
  --vcf chr22.vcf.gz \
  --out chr22_CEU.hmm.json \
  --population CEU
```

### Generate synthetic sequences

```bash
phantom-science generate \
  --hmm-model-path chr22_CEU.hmm.json \
  --population CEU \
  --reads 10000 \
  --output synthetic_CEU.fastq
```

### Measure indistinguishability (planned subcommand)

A `validate-indistinguishability` subcommand is planned for a future release.
It will run the Multinomial Naive Bayes cross-validation comparing real vs synthetic
sequences and report the measured accuracy.

---

## 7. References

1. **Weir & Cockerham (1984):** Estimating F-statistics for the analysis of population
   structure. *Evolution* 38(6), 1358–1370.

2. **Lewontin (1964):** The interaction of selection and linkage. *Genetics* 49(1), 49–67.

3. **Hill & Robertson (1968):** Linkage disequilibrium in finite populations.
   *Theoretical and Applied Genetics* 38(6), 226–231.

4. **Bird (1980):** DNA methylation and the frequency of CpG in animal DNA.
   *Nucleic Acids Research* 8(7), 1499–1504.

5. **Lander et al. (2001):** Initial sequencing and analysis of the human genome.
   *Nature* 409, 860–921.

6. **Auton et al. / 1000 Genomes Project Consortium (2015):** A global reference for
   human genetic variation. *Nature* 526, 68–74.

7. **Kelleher, Etheridge & McVean (2016):** Efficient coalescent simulation and
   genealogical analysis for large sample sizes. *PLOS Computational Biology* 12(5).

8. **Wigginton, Cutler & Abecasis (2005):** A note on exact tests of Hardy-Weinberg
   equilibrium. *Am. J. Hum. Genet.* 76, 887–893.

9. **McCallum & Nigam (1998):** A comparison of event models for Naive Bayes text
   classification. *AAAI Workshop on Learning for Text Categorization.*
