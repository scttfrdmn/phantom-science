// Package detmath provides the elementary functions and the normal and
// exponential random variates used on the output path, computed the same way
// on every CPU (docs/CHARTER.md, principle 4).
//
// The standard library does not guarantee this. math.Exp is assembly on
// amd64, arm64, loong64, riscv64 and s390x, and on amd64 it takes a fused
// multiply-add path only when the CPU has FMA. math.Log is assembly on amd64
// and s390x. Elsewhere these functions, and Pow, Log1p, Expm1, Lgamma, Erfc,
// Log2, Log10, Sin and Cos, are pure Go, which the compiler fuses into
// multiply-adds on arm64, ppc64le, riscv64 and s390x but not on default amd64.
// rand.(*Rand).NormFloat64 and ExpFloat64 call Exp and Log and contain fusable
// float32 arithmetic. Each of these can differ in the last bit between
// machines, and a last-bit difference changes output bytes whenever it crosses
// a rounding or rejection boundary.
//
// The functions here are the standard library's pure-Go algorithms (Go
// 1.27.2, from FreeBSD msun and Marsaglia and Tsang's ziggurat), ported
// without assembly and with every floating-point product wrapped in an
// explicit conversion so that it is rounded before it is added:
//
//	float64(x*y) + z // never fused
//
// Each result is therefore the IEEE 754 evaluation of the source expression,
// step by step, on every architecture and Go version. The results agree with
// the math package to within one ulp (see the tests). NormFloat64 and
// ExpFloat64 consume exactly the same draws from the generator as the rand
// methods, in the same order.
//
// The exact math functions (Sqrt, Abs, Floor, Ceil, Trunc, Round, Min, Max,
// Mod, Modf, Frexp, Ldexp, Copysign, the bit conversions and the special
// values) are correctly rounded or exact by definition and are used from the
// math package directly.
//
// internal/fmacheck checks this package on two levels: a source scan for
// products that could be fused, and a disassembly of the package compiled
// for arm64, amd64 (GOAMD64=v3), ppc64le, s390x, riscv64 and loong64 that
// fails on any fused multiply-add instruction. It also fails for any
// output-path use of the CPU-dependent math functions or rand methods.
//
// The ported files keep the Go authors' copyright notice; they are under the
// BSD-style license in LICENSE-GO in this directory (see also NOTICE).
package detmath
