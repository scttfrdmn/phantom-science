// Copyright 2009 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE-GO file in this directory.
//
// Ported to package detmath from Go 1.27.2 src/math/bits.go
// (the IEEE 754 binary64 layout constants only).

package detmath

const (
	mask     = 0x7FF
	shift    = 64 - 11 - 1
	bias     = 1023
	fracMask = 1<<shift - 1
)
