# Phantom Science Container Image

`docker/Dockerfile` builds a Phantom Science image with steganographic
fingerprinting. It needs no GPU, CUDA or credentials: the stegano library is
fetched from its public repository at the version pinned in `go.mod` and built
during the image build.

## Building

From the repository root:

```bash
docker build -f deploy/production/docker/Dockerfile -t phantom-science .
```

The build has two Debian bookworm stages:

1. **build** (`golang:1.26-bookworm`): installs Rust, runs
   `scripts/build-stegano.sh` to build `libstegano.so`, and builds the
   `phantom-science` binary with CGO linked against it.
2. **runtime** (`debian:bookworm-slim`): the binary in `/usr/local/bin`,
   `libstegano.so` in `/usr/local/lib`, and the example configurations from
   `configs/` in `/usr/share/phantom-science/configs`.

The container runs as an unprivileged user (UID/GID 65532) with `/data` as its
working directory. Its entrypoint is `phantom-science`, and with no arguments it
prints the help text. The health check runs `phantom-science --version`.

```bash
docker run --rm phantom-science --version
```

## Running

The container can only read and write the directories you mount, so mount a
host directory on `/data`. It must be writable by UID 65532, or you can run the
container with `--user "$(id -u):$(id -g)"`. Relative paths in a config, such as
`output_dir`, resolve against `/data`.

### Generating reads

Start from one of the example configs in the image (or in `configs/` in the
repository) and edit it for your run, for example `total_batches`,
`reads_per_batch` and `batch_interval`:

```bash
mkdir -p work && cd work
docker run --rm --entrypoint cat phantom-science \
  /usr/share/phantom-science/configs/sequencer/simple-run.json > my-run.json
# edit my-run.json

docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" \
  phantom-science generate --config my-run.json
```

The output lands in `work/` under the config's `output_dir`. A trained HMM
model can be placed in the same directory and passed with `--hmm-model-path`,
or a directory of models with `--model-dir` and `--population`.

### Fingerprinting and verifying

When the config enables the `quality_scores` fingerprinting layer, generation
needs a key. Either mount a key file (mode `0600`, named by `key_file` in the
config, relative to `/data`), or pass the key in the `PHANTOM_FINGERPRINT_KEY`
environment variable. Use the same key to verify. From `work/`:

```bash
openssl rand -base64 32 > fingerprint.key && chmod 600 fingerprint.key

docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" \
  -e PHANTOM_FINGERPRINT_KEY="$(cat fingerprint.key)" \
  phantom-science generate --config my-run.json

docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" \
  -e PHANTOM_FINGERPRINT_KEY="$(cat fingerprint.key)" \
  phantom-science verify sequencer_output/*.fastq
```

The glob in the `verify` example expands on the host, which works because the
host directory and `/data` have the same layout. `verify` also accepts
`--key-file`, `--customer` and `--dataset`; see `phantom-science verify --help`.

### Other commands

The image contains the same CLI as a local build: `generate`, `validate`,
`verify`, `train`, `scenario`, `transcriptome`, `signal`, `bench` and `tui`.
Run `phantom-science <command> --help` for each command's flags. `tui` needs an
interactive terminal (`docker run -it`). `signal` needs an ONT pore model
(`scripts/fetch-pore-model.sh`, run on the host) in the mounted directory.
Trio scenarios need a Python environment with msprime, which the image does
not have.
