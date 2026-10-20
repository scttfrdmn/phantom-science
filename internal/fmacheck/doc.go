// Package fmacheck holds a test that keeps generated output byte-identical
// across CPU architectures (docs/CHARTER.md, principle 4).
//
// The Go specification lets an implementation fuse x*y + z into a single
// fused multiply-add, which skips the rounding of x*y. The gc compiler does
// this on arm64, ppc64le, s390x and riscv64 (and on amd64 when built with
// GOAMD64=v3 or later), but not on default amd64, so the same seed would give
// different bytes on different machines. An explicit conversion rounds to the
// target precision and forbids the fusion:
//
//	r = float64(x*y) + z // never fused
//
// The test in this package type-checks the non-test source of every package
// whose floating-point arithmetic reaches generated output and fails, listing
// file:line, for each floating-point product that could be fused into an add
// or subtract. The package has no non-test code.
package fmacheck
