// Copyright 2009 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE-GO file in this directory.
//
// Ported to package detmath from Go 1.27.2 src/math/log10.go
// with every floating-point product explicitly rounded; see doc.go.

package detmath

import "math"

// Log10 returns the decimal logarithm of x.
// The special cases are the same as for [Log].
func Log10(x float64) float64 {
	return log10(x)
}

func log10(x float64) float64 {
	return float64(Log(x) * (1 / math.Ln10))
}

// Log2 returns the binary logarithm of x.
// The special cases are the same as for [Log].
func Log2(x float64) float64 {
	return log2(x)
}

func log2(x float64) float64 {
	frac, exp := math.Frexp(x)
	// Make sure exact powers of two give an exact answer.
	// Don't depend on Log(0.5)*(1/Ln2)+exp being exactly exp-1.
	if frac == 0.5 {
		return float64(exp - 1)
	}
	return float64(Log(frac)*(1/math.Ln2)) + float64(exp)
}
