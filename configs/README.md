# Configuration files

Two kinds of configuration are passed to `phantom-science`:

- a **sequencer config** (`--config`), which says what to generate: sample name,
  read length, paired or single end, reads per file, number of files, output
  directory, and so on;
- an optional **phantom config** (`-p` for `generate`, `validate` and the root
  command; `--phantom-config` for `tui`), which turns on fingerprinting and
  validation of the generated files.

Every sequencer config listed here was run with
`phantom-science generate --config FILE --total-reads 20 --seed 1 --output DIR`
and loads and generates, except the three v2 templates marked "description
only", which fail with an explicit error.

## Sequencer configs

The loader (`internal/cli/config_loader.go`) recognises four layouts:

| Layout | How it is recognised | Example |
|---|---|---|
| JSON | `.json` extension; keys are the JSON field names of `types.SequencerConfig` (`pkg/types/config.go`) | `configs/sequencer/simple-run.json` |
| Flat YAML | `.yaml`/`.yml` with the same keys as JSON | `examples/pacbio_revio_example.yaml` |
| v1 YAML | top-level `schema_version` | `configs/sequencer/v1.0/*.yaml` |
| v2 platform template | top-level `version: "2.x"` and a `specifications` section | `configs/sequencers/*.yaml` |

A file that sets none of the recognised keys is an error rather than an empty
config.

### Files

| File | Layout | Notes |
|---|---|---|
| `sequencer/simple-run.json` | JSON | 2×150, 10 files of 1M read pairs |
| `sequencer/miseq-v3.json`, `nextseq-2000-p3.json`, `novaseq-6000-s2.json` | JSON | instrument-named runs; see `sequencer/README.md` |
| `sequencer/v1.0/miseq-v3.yaml`, `nextseq-2000-p3.yaml`, `novaseq-6000-s2.yaml` | v1 YAML | the same runs in the v1 layout |
| `sequencers/ont_promethion24.yaml` | v2 template | Oxford Nanopore output (`platform: ont`): lognormal read lengths, SUP error model, MinKNOW headers, `sequencing_summary.txt`, gzip |
| `sequencers/pacbio_revio.yaml`, `element_aviti.yaml`, `mgi_dnbseq_t20x2.yaml` | v2 template | description only: loading fails with "only Oxford Nanopore templates are supported" |
| `dev/config.json` | JSON | 3 files of 10,000 read pairs with plaintext quality-score fingerprints (needs a stegano build) |
| `dev/fast-test.json` | JSON | 50 bp single end, 10,000 files of 100 reads |
| `dev/high-volume-test.json` | JSON | 1,000 files of 1,000 read pairs |
| `dev/million-files.json` | JSON | 500,000 files of 100 read pairs (one million files) |
| `dev/illumina-miseq.json`, `illumina-novaseq.json`, `nextseq-2000.json` | JSON | one large file pair each (2.5M, 50M and 10M read pairs) |

Use `--total-reads`, `--total-bases` or `--output` to override the size or
location without editing a file.

### Which fields have an effect

`generate` reads: `sample_name`, `lane`, `output_dir`, `read_length`,
`reads_per_batch` (reads or read pairs per file), `total_batches` (number of
files), `paired_end`, `base_quality`, `quality_decay`, `dataset_id`,
`run_id`, `flowcell_id`, `max_output_gb`/`max_output_tb` (a cap on the bases
of each flow cell), `enterprise.customer_id`, `fingerprinting` or
`enterprise.fingerprinting`, `hmm_model_path`, `population`, `model_dir`,
`platform`, `seed`, `read_length_distribution`, `error_model`, `ont`,
`output`, `total_bases`, `total_reads`, `flow_cells` and `workers`.

`sequencer_model` and `flow_cell_type` are used only by `tui` (run timing).
`batch_interval`, `q30_percentage`, `chemistry_version`, `instrument_id`,
`cluster_density_k_per_mm2`, `loading_concentration_pm`,
`phix_spike_in_percent`, `run_time_hours`, `max_cycles` and the
`*_channel_chemistry` flags are accepted but not read by the generator; they
describe the instrument. The vendor figures in these files have not been
checked against current specifications.

The v1 YAML layout maps only the older fields (see the comment at the top of
each `v1.0` file). For `seed`, `platform`, `read_length_distribution`,
`output` and the other newer settings, use JSON or flat YAML, or the
matching `generate` flags.

## Phantom configs

| File | Notes |
|---|---|
| `phantom/v1.0/enterprise-settings.yaml` | encrypted quality-score fingerprints, read-ID customer tags, and the built-in validation checks |

A phantom config is read for:

- `enterprise.customer_id` and `fingerprinting` (`enable`,
  `customer_fingerprint`, `dataset_fingerprint`, `steganography_layers`,
  `key_file`, `plaintext`). `quality_scores` is the only steganography layer
  that embeds a fingerprint; it needs the stegano library and an encryption
  key (from `key_file`, mode 600, or `PHANTOM_FINGERPRINT_KEY`) unless
  `plaintext: true` is set.
- `validation` (`enable` and the three test lists). Valid test names are
  `format_compliance`, `quality_distribution`, `sequence_complexity`,
  `read_length_consistency`, `hardy_weinberg`, `base_composition`,
  `gc_content`, `fastq_parsing`, `quality_encoding` and `header_format`;
  unknown names fail. With `enable: true` and no lists, the built-in suite of
  nine checks runs.

Other v1 keys (`audit_logging`, `regulatory_compliance`,
`academic_validation`, `security`, `cryptographic_signing`,
`tamper_detection`, `service_tier`, `compliance_level`, ...) are parsed but
have no effect, and phantom-science makes no compliance claims based on them.

## Usage

```bash
# Generate from a sequencer config
phantom-science generate --config configs/sequencer/v1.0/miseq-v3.yaml --total-reads 1000 --seed 1

# Oxford Nanopore reads
phantom-science generate --config configs/sequencers/ont_promethion24.yaml --total-bases 5M --seed 1

# With fingerprinting and validation (needs a stegano build and a key)
export PHANTOM_FINGERPRINT_KEY="$(openssl rand -base64 32)"
phantom-science generate --config configs/sequencer/simple-run.json \
  -p configs/phantom/v1.0/enterprise-settings.yaml --total-reads 1000 --seed 1
```

`internal/cli/config_loader_v2_test.go` checks that every
`configs/sequencers/*.yaml` and `configs/sequencer/v1.0/*.yaml` file either
loads or fails with an explicit error.
