# Project charter

These principles decide design questions in phantom-science. A change that
conflicts with one of them needs to change the principle first, in this file,
with a reason.

## 1. The data is the product

Generated reads, qualities, signal and metadata must be statistically
realistic for the platform they model, because their value is how they behave
in real tools: aligners, variant callers, basecallers, QC and classifiers.
Realism is measured against real data from public sources, with the
comparison recorded in the repository (see `docs/SCIENTIFIC_VALIDATION.md`
and `scripts/acceptance/`). Where output is not realistic, the documentation
says so.

## 2. Fingerprinting never alters the data

A fingerprint must not change any base, quality, length, signal sample or
other value of the generated data, and must not change their statistics.
Output with fingerprinting on must be indistinguishable from output with it
off, by any measurement of the data. Provenance is carried by how the data is
generated (a keyed seed), not by editing the data after generation. A
fingerprint that a user can see in the data is a defect: it degrades realism,
and it can be stripped.

Fingerprints are verified with the key and the run's records (manifest or
candidate IDs), not by extracting a payload from the file.

## 3. Claims are measured, not asserted

Every number in the documentation comes from a command or test in the
repository, or from a cited source. Code does not report fixed constants as
measurements. Unimplemented features are not described as working.

## 4. Output is reproducible

The same configuration and seed produce byte-identical output, whatever the
number of workers or the CPU architecture. Randomness comes only from seeded
streams, and floating-point products on the output path are explicitly
rounded so the compiler cannot fuse them (`internal/fmacheck` enforces this).
Known gap: the standard library's `math.Exp`/`Log` and `rand.NormFloat64` can
differ in the last bit between CPUs (#33).

## 5. Limitations are stated

What a model does not do, and where it has not been validated, is written
down next to what it does.
