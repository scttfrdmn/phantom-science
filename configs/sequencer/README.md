# Illumina sequencer configs

JSON sequencer configs named after Illumina instruments, plus their v1 YAML
equivalents in `v1.0/`. See `../README.md` for the layouts and for which
fields affect generation.

| File | Read length | Files × read pairs per file | Notes |
|---|---|---|---|
| `simple-run.json` | 2×150 | 10 × 1,000,000 | generic run for testing |
| `novaseq-6000-s2.json`, `v1.0/novaseq-6000-s2.yaml` | 2×150 | 4 × 20,000,000 | `max_output_tb: 3.0` caps the bases per flow cell |
| `nextseq-2000-p3.json`, `v1.0/nextseq-2000-p3.yaml` | 2×150 | 8 × 5,000,000 | `max_output_gb: 120` |
| `miseq-v3.json`, `v1.0/miseq-v3.yaml` | 2×300 | 6 × 2,500,000 | `max_output_gb: 15` |

The instrument name changes the output only through these values: all of
them produce the same Illumina-style FASTQ (headers
`@<sample>:<lane>:HXXXXXX:...`, files `<sample>_S1_L<lane>_R<mate>_<n>.fastq`).
Base qualities follow `base_quality` and `quality_decay`; there is no
per-instrument error or quality model. Fields such as `run_time_hours`,
`q30_percentage`, `cluster_density_k_per_mm2`, `loading_concentration_pm`,
`phix_spike_in_percent` and the chemistry flags describe the instrument and
are not read by the generator, and their values have not been checked against
current Illumina specifications.

## Usage

```bash
phantom-science generate --config configs/sequencer/novaseq-6000-s2.json --total-reads 10000 --seed 1
phantom-science generate --config configs/sequencer/v1.0/miseq-v3.yaml --total-reads 10000 --seed 1
phantom-science tui --config configs/sequencer/v1.0/nextseq-2000-p3.yaml
```

Without `--total-reads` or `--total-bases`, a run writes every file in the
config (for example 4 × 20M read pairs for `novaseq-6000-s2.json`).
