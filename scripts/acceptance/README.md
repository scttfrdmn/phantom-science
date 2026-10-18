# Acceptance scripts

These scripts check phantom-science output against the real tools of the
nanopore ecosystem. They need external software, so they are run by hand,
not in CI. Each one checks for its tools before doing any work and says what
is missing.

| Script | Checks | Issue | Status |
|---|---|---|---|
| `slow5_validate.sh` | BLOW5/SLOW5 files pass slow5tools and convert to pod5 | #15 | Passed: 132/132 checks with slow5tools 1.5.0 and blue-crab 0.6.0 (recorded in #15) |
| `dorado_identity.sh` | Dorado basecalls of synthetic signal match the source sequences | #14 | Not run yet |
| `transcriptome_bambu.sh` | Bambu quantification and novel-isoform discovery on simulated cDNA | #13 | Not run yet |
| `pycoqc.sh` | PycoQC reads a generated `sequencing_summary.txt` | #16 | Not run yet |
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
real data, e.g. a slice of `s3://ont-open-data/UHRR_HG002_2026.06`.

- Tools: `blue-crab`, `dorado`, `minimap2`, `samtools`, `awk` (required).
  `slow5tools` is optional. The level comparison also needs python3 with
  `pod5` and `numpy`.
- Arguments: `--blow5 FILE` (one file) and `--truth FASTA|FASTQ` (required).
  Optional: `--model` (default `sup`), `--real-pod5`, `--outdir` (default
  `./dorado_identity_out`), `--threads` (default 8), `--device`, `--rna`.
- Pass: median identity in [95%, 99.9%). Expect about 98-99% for R10.4.1
  with SUP. Identity of 99.9% or more means the noise model is too clean.
  Below 95% means the dwell, stall or skip parameters are off. Exits 0 on
  pass, 1 on fail, 2 on usage or tool errors.
- `signal` writes one BLOW5 file per `--reads-per-file` reads (default 4000).
  Either pass one chunk, raise `--reads-per-file` so the run fits in one
  file, or merge the chunks first:
  `slow5tools merge sig/ -o sig.blow5`.

## transcriptome_bambu.sh (#13)

Simulates ONT cDNA reads with `phantom-science transcriptome` from a genome
and annotation. It runs the standard long-read pipeline on them: Pychopper,
then `minimap2 -ax splice -uf -k14`, then Bambu. It then scores Bambu
against the simulation's truth (`counts.tsv`, `novel_isoforms.tsv`).

- Tools: `pychopper` >= 2.7, `minimap2` >= 2.22, `samtools` >= 1.10,
  `Rscript` with Bioconductor `bambu` >= 3.0. It also needs `go` to build
  phantom-science, unless `PHANTOM` points at a binary.
- Environment (defaults): `READS` (200000), `SEED` (1), `NOVEL_RATE` (0.05),
  `TIER` (sup), `KIT` (PCS114), `THREADS` (8), `NDR` (Bambu's default),
  `MIN_NOVEL_READS` (5), `PHANTOM`.
- Output: `OUTDIR/acceptance.tsv` and stdout. It reports the Pearson
  correlation of log1p(true reads) against log1p(Bambu counts) and the
  Spearman correlation, both over annotated transcripts. It also reports
  novel-isoform precision and recall, scored on injected exon-skip isoforms
  with at least `MIN_NOVEL_READS` reads. Alternative-3'-end isoforms keep
  their parent's intron chain, so they are counted but not scored.
- Pass criteria: the script exits 0 once the pipeline completes. Judge the
  metrics yourself; there is no hard threshold. Expect high count
  correlations, and recall that rises with the reads per novel isoform.

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
scripts/acceptance/dorado_identity.sh --blow5 sig/*_0.blow5 --truth sig/truth.fa --device metal
```

`truth.fa` holds, for each BLOW5 read ID, the sequence the signal was made
from. By default that is the read as written in the FASTQ, including the
simulated basecaller errors. If you have the error-free sequences, pass them
with `signal --truth true.fa` (keyed by read ID). The signal is then made
from the true molecule, and the identity measures Dorado alone. Non-UUID
read IDs are mapped to version-5 UUIDs; the original ID is kept in `truth.fa`
as `source=<id>`.
