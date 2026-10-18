// Package slow5 reads and writes nanopore raw-signal files in the SLOW5
// format: SLOW5 ASCII (.slow5) and its binary form BLOW5 (.blow5).
//
// The implementation follows the SLOW5 specification v0.2.0
// (Gamaarachchi et al., "SLOW5 Specification (version 0.2.0)", 14 May 2022,
// https://hasindu2008.github.io/slow5specs, docs/slow5-v0.2.0.pdf), cited
// below by section and table name, e.g. "spec, BLOW5 Header, Table 6". Where
// the specification is silent or ambiguous the behaviour of the reference
// implementation slow5lib (https://github.com/hasindu2008/slow5lib, MIT,
// src/slow5.c and src/slow5_press.c, library version 1.5.0) is followed and
// noted in a comment starting "slow5lib:".
//
// What is supported:
//
//   - SLOW5 ASCII read and write.
//   - BLOW5 read and write with record compression none, zlib or zstd
//     (spec codes 0, 1, 2) and signal compression none or svb-zd (codes 0, 1).
//     zstd uses github.com/klauspost/compress; svb-zd (zig-zag delta plus
//     StreamVByte) is implemented here in pure Go, byte-compatible with
//     slow5lib.
//   - All auxiliary field types of the specification: the eight signed and
//     unsigned integers, float, double, char, enum, char* and the 1D arrays.
//   - A chunking writer that spreads a run over several BLOW5 files named
//     <flowcell>_<runid>_<n>.blow5, mirroring PromethION's per-file read
//     counts.
//
// What is not: the slow5lib extension signal codec ex-zd (code 2, not part
// of specification v0.2.0) and the .idx index files (slow5tools index builds
// those).
//
// The package imports only the standard library and klauspost/compress, so
// any part of the genomics engine can use it without an import cycle.
package slow5
