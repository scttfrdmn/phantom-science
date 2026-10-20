// Package fmacheck holds the tests that keep generated output byte-identical
// across CPU architectures (docs/CHARTER.md, principle 4).
//
// The Go specification lets an implementation fuse x*y + z into a single
// fused multiply-add, which skips the rounding of x*y. The gc compiler does
// this on arm64, ppc64le, s390x, riscv64 and loong64 (and on amd64 when built
// with GOAMD64=v3 or later), but not on default amd64, so the same seed would
// give different bytes on different machines. An explicit conversion rounds
// to the target precision and forbids the fusion:
//
//	r = float64(x*y) + z // never fused
//
// The tests type-check the non-test source of internal/detmath and of every
// package whose floating-point arithmetic reaches generated output, and fail,
// listing file:line, for:
//
//   - each floating-point product that could be fused into an add or
//     subtract (TestNoFusableFloatProducts). A quotient by a power of two
//     counts as a product, since the compiler turns x/2 into x*0.5, and so
//     does a product converted to an unsigned integer, which some
//     architectures lower to a subtract;
//   - each use of a math function that is not exact on every CPU (Exp, Log,
//     Pow, Erfc, Lgamma, ...) or of rand's NormFloat64 and ExpFloat64
//     (TestNoCPUSpecificMath). The standard library computes these with
//     per-architecture assembly or fusable Go; use internal/detmath;
//   - each fused multiply-add instruction in the object code of those
//     packages compiled for arm64, amd64 (GOAMD64=v3), ppc64le, s390x,
//     riscv64 and loong64 (TestNoFusedInstructions), which checks the
//     compiler's actual output rather than the source.
//
// The package has no non-test code.
package fmacheck
