# syntax=docker/dockerfile:1
#
# Phantom Science container image with steganographic fingerprinting.
#
# Build from the repository root:
#   docker build -f deploy/production/docker/Dockerfile -t phantom-science .
#
# The build stage fetches and builds the stegano shared library at the version
# pinned in go.mod (scripts/build-stegano.sh), then links the Go binary against
# it with CGO. Both stages are Debian bookworm (glibc), so the library built in
# the first stage runs unchanged in the second. No GPU or CUDA is needed.

############################
# Build stage
############################
FROM golang:1.24-bookworm AS build

# The golang image already has git and a C toolchain; add Rust for stegano.
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl git \
    && rm -rf /var/lib/apt/lists/*
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
      | sh -s -- -y --no-modify-path --profile minimal --default-toolchain stable \
    && cargo --version

WORKDIR /src

# Download Go modules in their own layer so source edits do not refetch them.
COPY go.mod go.sum ./
RUN go mod download

COPY . .

# Build libstegano into /src/.stegano/target/release.
RUN scripts/build-stegano.sh

# Link against the freshly built library. At runtime it is found in
# /usr/local/lib (rpath, plus ldconfig in the runtime stage).
RUN CGO_ENABLED=1 \
    CGO_LDFLAGS="-L/src/.stegano/target/release -Wl,-rpath,/usr/local/lib" \
    go build -trimpath -o /out/bin/phantom-science . \
    && mkdir -p /out/lib \
    && cp /src/.stegano/target/release/libstegano.so /out/lib/

############################
# Runtime stage
############################
FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /out/lib/libstegano.so /usr/local/lib/libstegano.so
RUN ldconfig

COPY --from=build /out/bin/phantom-science /usr/local/bin/phantom-science

# Example sequencer and phantom configurations, for use with --config / -p.
COPY configs/ /usr/share/phantom-science/configs/

# Unprivileged user. /data is the working directory and the place to mount
# configs, keys and output directories.
RUN groupadd --gid 65532 phantom \
    && useradd --uid 65532 --gid 65532 --home-dir /data --no-create-home \
         --shell /usr/sbin/nologin phantom \
    && mkdir -p /data \
    && chown 65532:65532 /data
USER 65532:65532
WORKDIR /data

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD ["phantom-science", "--version"]

ENTRYPOINT ["phantom-science"]
CMD ["--help"]
