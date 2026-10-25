# Acceptance scripts

These scripts check phantom-science output against the real tools of the
nanopore ecosystem. They need external software, so they are run by hand,
not in CI. Each one checks for its tools before doing any work and says what
is missing.

| Script | Checks | Issue | Status |
|---|---|---|---|
| `slow5_validate.sh` | BLOW5/SLOW5 files pass slow5tools and convert to pod5 | #15 | Passed: 132/132 checks with slow5tools 1.5.0 and blue-crab 0.6.0 (recorded in #15) |
| `dorado_identity.sh` | Dorado basecalls of synthetic signal match the source sequences | #14 | Passed: 90.79% median identity against a squigulator baseline of 91.47% (threshold 90.46%), Dorado 2.1.2 SUP v5.2.0 |
| `calibrate/` | Ablation of the synthetic-vs-real basecalling gap, DNA R10.4.1 and RNA004 | #28 | Run 2026-10-20: results in `calibrate/results/` (see below) |
| `transcriptome_bambu.sh` | Bambu quantification and novel-isoform discovery on simulated cDNA | #13 | Completed on GRCh38 chr22 / GENCODE v49: count Pearson 0.872, Spearman 0.832; novel precision 0.90, recall 0.55 |
| `pycoqc.sh` | PycoQC reads a generated `sequencing_summary.txt` | #16 | Passed: PycoQC 2.5.2 counted all 537 reads |
| `illumina_heldout.sh` | The Illumina run profiles (three single-run, one pooled) against three held-out NovaSeq 6000 runs and against real-run-vs-real-run distances | #35, #36 | Run 2026-10-25: results in `calibrate/illumina/results/heldout/` (see below) |
| `illumina_indels.sh` | Illumina indel error rates where mates overlap (three fitted runs, their pool, three held-out runs), the fit of the indel model, and the closed loop on generated reads | #30 | Run 2026-10-25: results in `calibrate/illumina/results/indels/` (see below) |
| `../fetch-pore-model.sh` | Downloads an ONT k-mer pore model for `signal` | #14 | Used for the end-to-end example below |

Until a script has been run, the feature it checks has only the repository's
own tests behind it. S3 output (#16) is tested only against a fake uploader,
not against real AWS.

Commands below are run from the repository root. Build the binary first with
`make build` or `CGO_ENABLED=0 go build -o phantom-science .`.

## fetch-pore-model.sh

Downloads an ONT k-mer level table from
[nanoporetech/kmer_models](https://github.com/nanoporetech/kmer_models). The
tables are not shipped with phantom-science; read ONT's licence in that
repository before you use or share them.

- Tools: `curl` or `wget`.
- Models: `dna_r10.4.1_e8.2_400bps` (default), `dna_r10.4.1_e8.2_260bps`, `rna004`
  (`--list` prints them).
- Output: `DEST_DIR/<model>/9mer_levels_v1.txt`. The path is printed on
  stdout. `DEST_DIR` defaults to
  `${PHANTOM_PORE_MODEL_DIR:-$HOME/.cache/phantom-science/pore-models}`.
  `KMER_MODELS_REF` selects a branch, tag or commit (default `master`).
- Pass: the file has exactly 262144 well-formed 9-mer rows. Exits 1 if the
  download fails or the table is malformed, 2 for an unknown model or no
  downloader.

```sh
model=$(scripts/fetch-pore-model.sh dna_r10.4.1_e8.2_400bps)
```

These are two-column tables in standardised units. `phantom-science signal`
detects that and scales them to pA with the profile's scaling.

## slow5_validate.sh (#15)

Validates BLOW5/SLOW5 files with the reference tools. For each file it runs
`slow5tools quickcheck`, `view` (full decode), `stats`, `index` and `get`
(first and last read). It then re-encodes the file to zlib + svb-zd and
checks the records come out unchanged. Every directory holding more than
one `.blow5` file must also `merge` into one valid file. If blue-crab is
installed, each BLOW5 file also goes BLOW5 -> pod5 -> BLOW5, and that round
trip must keep every read ID and sample. Zstd files are re-encoded to zlib
first, because the PyPI blue-crab build cannot read zstd.

- Tools: `slow5tools` (required; build with `make zstd=1`,
  https://github.com/hasindu2008/slow5tools), `blue-crab` and `pod5`
  (optional; `pip install blue-crab`, Python 3.8-3.12). You can give paths
  with `SLOW5TOOLS`, `BLUECRAB` and `POD5`. slow5tools has no `s2p` command.
  BLOW5 -> pod5 conversion is blue-crab's.
- Without arguments it first writes sample files with
  `go test -run TestWriteSampleFiles ./internal/genomics/slow5/`, so it also
  needs `go`.
- Exit status: 0 all checks passed; 1 a check failed; 2 slow5tools missing;
  3 slow5tools checks passed but the pod5 step was skipped (no blue-crab).

```sh
scripts/acceptance/slow5_validate.sh                 # writer sample files
scripts/acceptance/slow5_validate.sh sig/            # output of `phantom-science signal`
SLOW5TOOLS=/opt/slow5tools/slow5tools BLUECRAB=$VENV/bin/blue-crab \
    scripts/acceptance/slow5_validate.sh sig/*.blow5
```

## dorado_identity.sh (#14)

Converts synthetic BLOW5 to pod5 (`blue-crab s2p`; if slow5tools is present
it first re-encodes to zlib + svb-zd). It then basecalls with Dorado
(`--emit-moves`), aligns the calls to the truth with minimap2 (`map-ont`, or
`splice -uf -k14` with `--rna`) and reports median, mean and aggregate
identity (PAF matches / alignment block length, primary alignments) and
mean query coverage. It also prints a dwell histogram from the move tables.
With `--real-pod5` it compares dwell and pA level distributions against
real data, e.g. a pod5 file from `s3://ont-open-data/pgx_as_2025.07`
(PromethION FLO-PRO114M, 5 kHz).

- Tools: `blue-crab`, `dorado`, `minimap2`, `samtools`, `awk` (required).
  `slow5tools` and `squigulator` are optional. The level comparison also
  needs python3 with `pod5` and `numpy`.
- Arguments: `--blow5 FILE` (one file) and `--truth FASTA|FASTQ` (required).
  Optional: `--model` (default `sup`), `--real-pod5`, `--outdir` (default
  `./dorado_identity_out`), `--threads` (default 8), `--device`, `--rna`,
  `--squigulator [PATH]`, `--tolerance` (default 1.0 points),
  `--min-identity` (default 90.5%).
- Pass: median identity at or above the threshold and below 99.9% (at or
  above 99.9% means the noise model is too clean). Exits 0 on pass, 1 on
  fail, 2 on usage or tool errors.
- Threshold: synthetic signal from ONT's k-mer level tables cannot reach the
  identity of real reads. On 2026-10-18, with Dorado 2.1.2 and
  `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`, real R10.4.1 reads
  (`pgx_as_2025.07`, 348 reads at MAPQ >= 20 against GRCh38) reached 98.52%
  median identity, and squigulator 0.5.0, the algorithm `signal`
  reimplements, reached 91.47% on the same truth as phantom-science. The
  pass threshold is therefore relative to squigulator: with
  `--squigulator`, the script simulates the same truth with squigulator,
  basecalls it the same way, and passes `signal` output within `--tolerance`
  points of it. Without it, `--min-identity` applies (90.5%: that baseline
  minus one point).
- Measured with the default `dna_r10.4.1_e8.2_400bps_5khz` profile, on the
  end-to-end example below (158 reads, seed 1): median identity 90.79%;
  dwell from Dorado's moves mean 13.03 samples/base, CV 0.886 (real: 13.12,
  CV 1.077 over 5,325 reads); signal mean 93.8 pA, SD 23.8 (real: 97.8 pA,
  SD 24.3).
- `signal --profile-set NAME=VALUE` overrides profile parameters for
  calibration runs; the overrides are recorded in the manifest.
- `signal` writes one BLOW5 file per `--reads-per-file` reads (default 4000).
  Either pass one chunk, raise `--reads-per-file` so the run fits in one
  file, or merge the chunks first:
  `slow5tools merge sig/ -o sig.blow5`.

## calibrate/ (#28)

Measures where the identity gap between synthetic and real signal comes
from. Each condition is a truth read set plus a simulator (`signal` with
`--profile-set` overrides, squigulator, or a real pod5). `eval.sh` converts
the signal to pod5, basecalls it with Dorado (`--emit-moves`), aligns the
calls to the truth with minimap2 and appends one row to a results TSV:
median identity with a bootstrap 95% CI over reads (1000 resamples), Dorado
move-table dwell statistics, signal statistics measured the same way on
synthetic and real reads (robust per-sample noise, lag-1 autocorrelation
of first differences, per-base level spread, and the correlation of
per-base levels with the k-mer table), the git SHA, the `TOOLS.lock`
sha256, the device and the date. Intermediates are deleted once a row is
written.

- Files: `conditions.tsv` (DNA) and `conditions_rna.tsv` (RNA004) list the
  conditions and pin the inputs by sha256; `run_conditions.sh` runs them in
  order and skips rows already in the results; `make_inputs.py` builds the
  truth sets (`regions.tsv` pins the 2000 GRCh38 regions); `measure.py`
  scores one condition; `squig_table.py` recovers squigulator's built-in
  level tables from `--ideal`/`--ideal-time` signal of a de Bruijn
  sequence; `fit_table.py` and `kmer_noise.py` compare per-k-mer levels and
  noise between real and synthetic reads.
- Tools: those of `dorado_identity.sh` plus squigulator, `pod5` and numpy
  (the blue-crab venv). Paths default to the `/tmp/ps-accept` layout; see
  the script headers.
- Inputs: DNA, 2000 GRCh38 pseudo-reads of 2-8 kb (10.0 Mb, GC 40.9%) and
  2000 uniform random reads of the same lengths from `generate`; real reads
  from `pgx_as_2025.07` (all 5325, MAPQ >= 20 against GRCh38). RNA004, 1200
  reads of the 6 control constructs of `rna-modbase-validation_2025.03`
  and the first 3000 reads of its `control_rep1.pod5`.

DNA, Dorado 2.1.2 `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` on Metal
(`results/dna.tsv`; the profile is the `dna_r10.4.1_e8.2_400bps_5khz`
preset unless stated):

| Condition | Median identity [95% CI] |
|---|---|
| Real reads | 98.54% [98.52, 98.58] |
| Uniform random reads | 90.80% [90.76, 90.82] |
| GRCh38 reads | 91.75% [91.67, 91.82] |
| squigulator 0.5.0, uniform / GRCh38 | 91.50% / 93.74% |
| squigulator with ONT's table, GRCh38 | 91.63% |
| GRCh38, squigulator's levels (constant stdv) | 93.81% [93.74, 93.87] |
| ... and no stalls | 93.96% [93.90, 94.00] |
| GRCh38, squigulator's levels and per-k-mer stdv | 92.88% |
| `level_tau` 0.5 / 1 / 2 / 4 samples | 91.49% / 91.29% / 87.77% / 72.90% |
| `level_tau` 8 | 254 of 2000 reads aligned |
| `level_noise_scale` 1.2 / 0.5 | 91.59% / 5 of 2000 aligned |
| +5.2 pA level shift (`model_shift` 98.921) | 91.71% |
| ONT's table refit to real reads (`fit_table.py`) | 77.54% |

- Composition is worth about 1 point with ONT's table (2 with
  squigulator's levels), and the level table about 2 points; everything
  else is worth at most 0.2 points or makes identity worse. Level smoothing
  lowers identity even at τ = 0.5 samples, where the autocorrelation of
  first differences already matches real reads (+0.046 vs +0.048).
- The remaining 4.6 points are shared with squigulator and are the k-mer
  level model: real per-base levels fit the tables at r = 0.70 (ONT) and
  0.72 (squigulator) against 0.85 for synthetic signal, and CpG k-mers sit
  3.8 pA lower in native, methylated DNA. Refitting the table from
  Dorado's move segmentation made things worse; that needs a resquiggler.
- Dwell from Dorado's moves: real 13.12 samples/base, CV 1.08; the preset
  gives 13.12, CV 0.94. Without stalls the CV drops to 0.61, so stalls stay.
- squigulator's built-in R10 table is a different table from ONT's
  `9mer_levels_v1` (r = 0.994, 2.6 pA RMS apart after a linear fit). Its
  RNA004 table is ONT's `rna004` table at 17.57 x + 84.04 pA with a
  per-k-mer stdv. Pass either to `signal --pore-model` after recovering it
  with `squig_table.py`.

RNA004, Dorado 2.1.2 `rna004_130bps_sup@v5.3.0` (`results/rna.tsv`;
aligned with `minimap2 -x map-ont -k5 -w1 -m10 -s15 -n2`, because at
`-k10 -w5` only about 1% of synthetic reads of the 100-nt constructs
align). Columns: aligned reads, median identity of aligned reads, and the
median over all reads of matched bases / construct length:

| Condition | Aligned | Median identity | Matched |
|---|---|---|---|
| Real reads | 99.9% | 97.67% [97.65, 97.70] | 0.84 |
| squigulator 0.5.0 | 97.9% | 82.86% | 0.39 |
| Old `rna004_130bps_4khz` preset | 61.7% | 63.41% | 0.12 |
| Calibrated preset (`rna_cal`) | 72.1% | 72.92% [72.09, 73.75] | 0.28 |
| ... with squigulator's per-k-mer stdv | 96.8% | 77.78% | 0.62 |

The calibrated values (adapter 0.7 s at 62 pA, Dorado's 17/80 pA scaling,
3.0 pA per-sample noise and no white noise or drift, dwell CV 0.31, no
skips) are measured on the real reads; `rna_preset` reruns the committed
preset without overrides and reproduces `rna_cal` exactly. Within the CIs,
speed (115 instead of 130 bps: 72.08%), no stalls (73.33%, but dwell CV
0.60 against 1.15 for real reads) and level smoothing (`level_tau` 0.5:
73.17%) do not change identity. The main measured lever is a per-k-mer
stdv, which ONT's RNA004 table does not carry.

```sh
scripts/acceptance/calibrate/run_conditions.sh --inputs          # DNA
scripts/acceptance/calibrate/run_conditions.sh \
    --conditions scripts/acceptance/calibrate/conditions_rna.tsv \
    --results scripts/acceptance/calibrate/results/rna.tsv       # RNA004
```

## illumina_heldout.sh (#35)

```sh
scripts/acceptance/illumina_heldout.sh hg38.primary.fa   # or no argument: downloads UCSC hg38
```

Streams the first 200,000 read pairs of each run in `RUNS` from ENA
(rejecting runs without binned `#,:F` qualities): the fitted runs DRR276626,
DRR394916 and ERR15600747 and the held-out runs DRR260182, SRR30646135 and
ERR12737258. For each profile in `PROFILES` (`NAME=RUN`, or
`NAME=RUN+RUN+...` for a pooled profile, whose runs all count as its own) it
runs `run_e2e_tests` paired against each run (real vs real, real vs
generated), and for every pair of runs `-real-only -real-fastq-b`, and
compares calibration statistics (`illumina_stats.py`, minimap2 alignments,
`heldout_compare.py`) of every run with the generated reads. Results, with
the md5 of the streamed reads, are in `calibrate/illumina/results/heldout/`:
`classifier_matrix.tsv`, `classifier_spread.tsv`, `profile_matrix.tsv`,
`profile_matching.tsv`, `profile_spread.tsv`, `calibration_stats.tsv`,
`cycle_delta.tsv`, `structure_compare.tsv`, `runs.tsv`, `run.tsv`, and per
read set `sets/NAME/`.
About 45 minutes on 12 cores; 9 GB for the minimap2 index.

## illumina_indels.sh (#30)

```sh
scripts/acceptance/illumina_indels.sh hg38.primary.fa
```

Measures indel errors where the two mates of a proper pair overlap: an
indel in one mate only, after left-normalisation against the reference
(rebuilt from SEQ, CIGAR and MD), is a sequencing error; one in both is the
sample's (`illumina_stats.py indels`). Left out as alignment artefacts:
mates that end inside the repeat, an indel of the other mate within 10
bases or the other mate's mismatches there (the same difference aligned
another way), an insertion and a deletion of one mate within 10 bases
(mismatches aligned as two gaps), and errors of different pairs within 50
bases (misaligned repeats). Errors are counted per homopolymer run by its
length, per short tandem repeat (2-6 base units) by its copies and per
base, with exact Poisson 95 % intervals. Only the overlap is measured (from
cycle insert − read length to the end of both mates, so early cycles come
only from short inserts), and indels within 10 bases of an alignment end
are not compared.

The runs of `RUNS` (default DRR394916 and DRR276626, 4,000,000 pairs each,
whose long inserts rarely overlap, and ERR15600747, 1,000,000) and
`HELDOUT` (DRR260182, SRR30646135, ERR12737258, 1,000,000 each; no profile
is fitted to them) are streamed from ENA and aligned with minimap2; the
pool is the first 1,000,000 pairs of each fitted run. `fit_indels.py` fits
the model's shape to the pool and each profile's insertion and deletion
rates to its run. With `RECOVERY` (default
`novaseq6000-err15600747=2000000`) it builds phantom-science with that
model, generates reads from the reference (seed `SEED` + 1), measures them
(`recovery-PROFILE`) and refits with the rates multiplied by the run's
measured rate over the generated reads' (×1.18 insertions, ×1.06
deletions). It writes `internal/genomics/illumina/profile_indels.go` and
`indel_fit.tsv` (the fits, the held-out runs against the pooled model, the
recovery, and tests of a cycle or mate dependence). Then the closed loop:
reads generated with each profile of `GENERATED` (default
`novaseq6000-pooled=3000000 novaseq6000-err15600747=2000000`) measured
the same way (`generated-PROFILE`, checked by
`TestCalibration_IndelClosedLoop`), and of `FLOOR` without indel errors
(`floor-PROFILE`: what the method finds where there is nothing). Results,
with the md5 of the streamed and generated reads, are in
`calibrate/illumina/results/indels/`: `summary.tsv` (every set's rates),
`indel_fit.tsv`, `runs.tsv`, `run.tsv` and per set `indel_*.tsv`
(summary, classes by repeat, coarser `indel_repeats.tsv`, lengths, sizes,
cycles, qualities, the sample's indels, every error in `indel_events.tsv`).
About 70 minutes on 12 cores with the downloads; 9 GB for the minimap2
index.

| Set | Pairs | Insertions per base (95 %), errors | Deletions per base (95 %), errors |
|---|---|---|---|
| pool (3 × 1,000,000) | 3,000,000 | 4.6·10⁻⁶ (3.8-5.5), 118 | 7.5·10⁻⁶ (6.5-8.7), 194 |
| DRR394916 | 4,000,000 | 6.5·10⁻⁶ (4.3-9.5), 28 | 9.1·10⁻⁶ (6.5-12.5), 39 |
| DRR276626 | 4,000,000 | 6.8·10⁻⁶ (5.0-9.1), 47 | 8.7·10⁻⁶ (6.7-11.3), 60 |
| ERR15600747 | 1,000,000 | 4.5·10⁻⁶ (3.6-5.4), 102 | 7.8·10⁻⁶ (6.7-9.0), 178 |
| DRR260182 (held out) | 1,000,000 | 4.8·10⁻⁶ (4.1-5.7), 145 | 4.4·10⁻⁶ (3.7-5.2), 132 |
| SRR30646135 (held out) | 1,000,000 | 5.4·10⁻⁶ (3.0-8.8), 15 | 7.5·10⁻⁶ (4.6-11.5), 21 |
| ERR12737258 (held out) | 1,000,000 | 4.6·10⁻⁶ (3.6-5.7), 79 | 5.2·10⁻⁶ (4.2-6.4), 90 |
| generated, novaseq6000-pooled | 3,000,000 | 4.6·10⁻⁶ (3.4-6.1), 50 | 7.4·10⁻⁶ (5.9-9.2), 81 |
| generated, novaseq6000-err15600747 | 2,000,000 | 5.2·10⁻⁶ (4.6-6.0), 230 | 7.8·10⁻⁶ (7.0-8.7), 342 |
| generated without indel errors (floor) | 2,000,000 | 1.1·10⁻⁷ (0.4-2.7), 5 | 2.0·10⁻⁷ (0.9-3.9), 9 |

## transcriptome_bambu.sh (#13)

Simulates ONT cDNA reads with `phantom-science transcriptome` from a genome
and annotation. It runs the standard long-read pipeline on them: Pychopper,
then `minimap2 -ax splice -uf -k14`, then Bambu. It then scores Bambu
against the simulation's truth (`counts.tsv`, `novel_isoforms.tsv`).

- Tools: `pychopper` >= 2.7, `minimap2` >= 2.22, `samtools` >= 1.10,
  `Rscript` with Bioconductor `bambu` >= 3.0. It also needs `go` to build
  phantom-science, unless `PHANTOM` points at a binary.
- Environment (defaults): `READS` (200000), `SEED` (1), `NOVEL_RATE` (0.2),
  `TIER` (sup), `KIT` (PCS114), `THREADS` (8), `NDR` (1), `MIN_NOVEL_READS`
  (5), `MIN_GENE_FRACTION` (0.05), `PHANTOM`.
- Output: `OUTDIR/acceptance.tsv` and stdout. It reports the Pearson
  correlation of log1p(true reads) against log1p(Bambu counts) and the
  Spearman correlation, both over annotated transcripts. It also reports
  novel-isoform precision and recall, scored on injected exon-skip isoforms.
  Alternative-3'-end isoforms keep their parent's intron chain, so they are
  counted but not scored.
- Detectable isoforms: recall is scored over injected exon-skip isoforms
  with at least `MIN_NOVEL_READS` reads and at least `MIN_GENE_FRACTION` of
  their gene's reads, because Bambu drops novel candidates below 5% of their
  gene's reads (`min.readFractionByGene`). `novel_recall_any_fraction` also
  reports recall over every isoform with enough reads. `NOVEL_RATE` 0.2
  gives injected isoforms a share Bambu can find; at 0.05 almost none reach
  5% of their gene.
- `NDR` 1 turns off Bambu's novel-discovery-rate filter, which is tuned for
  real data where most candidates are false. Here the injected isoforms are
  the only true novel transcripts, so precision still counts false
  discoveries. `NDR=` (empty) uses Bambu's default: on the run below that
  gave precision 0.90 and recall 0.135.
- Measured on 2026-10-18 (GRCh38 chr22 from UCSC, GENCODE v49 chr22, 11,615
  transcripts, defaults above; Pychopper 2.7.10 found primers in 99.86% of
  reads): count Pearson 0.872, Spearman 0.832; Bambu 3.8.3 reported 335
  novel multi-exon isoforms, precision 0.90 (302 are injected exon skips);
  recall 0.55 over 319 detectable isoforms, 0.23 over the 763 with at least
  5 reads. Of the 265 new junctions of the detectable isoforms, 228 have at
  least 5 aligned reads, so the misses are Bambu's read-class filters, not
  missing evidence.
- Pass criteria: the script exits 0 once the pipeline completes. There is
  no hard threshold; compare with the figures above.

```sh
awk '$1=="chr22"' gencode.v46.primary_assembly.annotation.gtf > chr22.gtf
scripts/acceptance/transcriptome_bambu.sh chr22.fa chr22.gtf bambu_out
```

## pycoqc.sh (#16)

Runs PycoQC on a `sequencing_summary.txt` written by `generate` for an ONT
run. The report must count every read of the summary.

- Tools: `pycoQC` (`pip install pycoQC`) and `python3`. Without arguments it
  also needs the phantom-science binary (`PHANTOM_BIN`, default
  `phantom-science` on PATH) and an ONT config (`PHANTOM_ONT_CONFIG`, default
  `configs/sequencers/ont_promethion24.yaml`). It then generates a 20 Mb run
  with `--seed 16 --workers 4` in a temporary directory.
- Arguments: `[summary-or-output-dir] [report-dir]`. Given a directory, it
  uses the first `sequencing_summary*.txt` under it.
- Pass: PycoQC exits 0, writes a JSON report, and its
  `All Reads.basecall.reads_number` equals the summary's row count. The HTML
  report is kept in `report-dir`.

```sh
PHANTOM_BIN=./phantom-science scripts/acceptance/pycoqc.sh
scripts/acceptance/pycoqc.sh reads/ qc_report/
```

## Tool versions

These combinations worked on macOS arm64 (2026-10-18). Newer library
versions broke each tool:

- PycoQC 2.5.2 needs `pandas==1.5.3`, `numpy<2`, `plotly<6` and
  `setuptools<81` (Python 3.11). With pandas 2 it fails with "sort array is
  read-only".
- Pychopper 2.7.10 needs `pandas<3`; its plotting step fails on pandas 3.
- Bambu 3.8.3 needs xgboost < 2 (xgboost 3 fails with "Got numeric 'y'").
  Bioconda has no osx-arm64 Bambu, so install it from Bioconductor into a
  conda env with R 4.4 and `compilers`, after building xgboost 1.7.11.1 from
  the CRAN archive.
- Dorado 2.1.2 (`dorado-2.1.2-osx-arm64.zip` from ONT's CDN), blue-crab
  0.6.0, squigulator 0.5.0 (built from source), minimap2 2.31.

## End to end: reads -> signal -> basecalling identity

```sh
# 1. ONT reads with MinKNOW headers (UUID read IDs, runid, ch, start_time ...).
./phantom-science generate --config configs/sequencers/ont_promethion24.yaml \
    --platform ont --total-bases 5M --seed 1 --output reads

# 2. A pore model (not shipped; ONT licence applies).
model=$(scripts/fetch-pore-model.sh dna_r10.4.1_e8.2_400bps)

# 3. Raw signal: BLOW5 (zlib + svb-zd), truth.fa and phantom-manifest.json.
#    Read IDs and run metadata are carried over from the FASTQ headers.
./phantom-science signal --input reads/*.fastq* --out sig --pore-model "$model" \
    --seed 1 --reads-per-file 100000

# 4. Validate the files, then basecall and score.
scripts/acceptance/slow5_validate.sh sig/
scripts/acceptance/dorado_identity.sh --blow5 sig/*_0.blow5 --truth sig/truth.fa --device metal \
    --squigulator
```

`truth.fa` holds, for each BLOW5 read ID, the sequence the signal was made
from. By default that is the read as written in the FASTQ, including the
simulated basecaller errors. If you have the error-free sequences, pass them
with `signal --truth true.fa` (keyed by read ID). The signal is then made
from the true molecule, and the identity measures Dorado alone. Non-UUID
read IDs are mapped to version-5 UUIDs; the original ID is kept in `truth.fa`
as `source=<id>`.
