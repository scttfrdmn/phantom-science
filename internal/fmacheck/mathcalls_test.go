//go:build !release

package fmacheck

import (
	"fmt"
	"go/ast"
	"go/importer"
	"go/parser"
	"go/token"
	"go/types"
	"sort"
	"strings"
	"testing"
)

// exactMath are the math functions whose results are exact or correctly
// rounded by definition, so every correct implementation (assembly or Go, on
// any CPU) returns the same bits. Everything else in package math is an
// approximation whose last bit depends on the implementation; the output
// path must use internal/detmath instead.
var exactMath = map[string]bool{
	"Abs": true, "Ceil": true, "Copysign": true, "Dim": true,
	"FMA": true, "Float32bits": true, "Float32frombits": true,
	"Float64bits": true, "Float64frombits": true, "Floor": true,
	"Frexp": true, "Inf": true, "IsInf": true, "IsNaN": true,
	"Ldexp": true, "Max": true, "Min": true, "Mod": true, "Modf": true,
	"NaN": true, "Nextafter": true, "Nextafter32": true, "Remainder": true,
	"Round": true, "RoundToEven": true, "Signbit": true, "Sqrt": true,
	"Trunc": true,
}

// cpuSpecificRand are the math/rand and math/rand/v2 functions and methods
// that compute their result with math.Exp or math.Log (and fusable float32
// arithmetic). Their integer and Float64 draws are exact and allowed.
var cpuSpecificRand = map[string]string{
	"NormFloat64": "detmath.NormFloat64(r)",
	"ExpFloat64":  "detmath.ExpFloat64(r)",
	"NewZipf":     "an inverse-CDF sampler built on detmath",
}

// detmathPackage is exempt from the call check: it is the implementation.
const detmathPackage = "/internal/detmath"

// TestNoCPUSpecificMath fails for every use, on the output path, of a math
// function that is not exact on every CPU, or of the rand methods built on
// them. Use the internal/detmath equivalent, which computes the same
// algorithm with every product explicitly rounded.
func TestNoCPUSpecificMath(t *testing.T) {
	pkgs := loadPackages(t, outputPackages)
	var all []string
	for _, p := range pkgs {
		if strings.HasSuffix(p.listed.ImportPath, detmathPackage) {
			continue
		}
		all = append(all, cpuSpecificUses(p)...)
	}
	for _, f := range all {
		t.Error(f)
	}
	if len(all) > 0 {
		t.Logf("%d CPU-specific math or rand uses on the output path; see internal/detmath", len(all))
	}
}

// cpuSpecificUses lists, as "file:line: message", every reference in p to a
// non-exact math function or a CPU-specific rand function or method,
// whether it is called or taken as a function value.
func cpuSpecificUses(p *loaded) []string {
	var out []string
	for id, obj := range p.info.Uses {
		fn, ok := obj.(*types.Func)
		if !ok || fn.Pkg() == nil {
			continue
		}
		pos := p.fset.Position(id.Pos())
		switch fn.Pkg().Path() {
		case "math":
			if !exactMath[fn.Name()] {
				out = append(out, fmt.Sprintf("%s: math.%s differs between CPUs; use detmath.%s",
					pos, fn.Name(), fn.Name()))
			}
		case "math/rand", "math/rand/v2":
			if repl, bad := cpuSpecificRand[fn.Name()]; bad {
				out = append(out, fmt.Sprintf("%s: %s.%s uses math.Exp/Log and differs between CPUs; use %s",
					pos, fn.Pkg().Path(), fn.Name(), repl))
			}
		case "math/cmplx":
			out = append(out, fmt.Sprintf("%s: cmplx.%s is built on math.Exp/Log/Sin; not portable",
				pos, fn.Name()))
		}
	}
	sort.Strings(out)
	return out
}

// TestCPUSpecificUsesFound runs the call check on a small source file.
func TestCPUSpecificUsesFound(t *testing.T) {
	const src = `package p

import (
	"math"
	"math/rand/v2"
	mr "math/rand"
)

func f(r *rand.Rand, r1 *mr.Rand, x float64) float64 {
	a := math.Exp(x)        // flag
	b := math.Sqrt(x)       // safe
	c := r.NormFloat64()    // flag
	d := r.Float64()        // safe
	e := r1.ExpFloat64()    // flag
	g := math.Floor(x)      // safe
	h := math.Log           // flag: function value
	l, _ := math.Lgamma(x)  // flag
	n := rand.NormFloat64() // flag: package function
	m := math.Pi * math.Abs(x)
	return a + b + c + d + e + g + h(x) + l + n + m
}
`
	fset := token.NewFileSet()
	file, err := parser.ParseFile(fset, "p.go", src, 0)
	if err != nil {
		t.Fatal(err)
	}
	info := &types.Info{
		Types: map[ast.Expr]types.TypeAndValue{},
		Defs:  map[*ast.Ident]types.Object{},
		Uses:  map[*ast.Ident]types.Object{},
	}
	conf := types.Config{Importer: importer.ForCompiler(fset, "source", nil)}
	pkg, err := conf.Check("example.com/p", fset, []*ast.File{file}, info)
	if err != nil {
		t.Fatal(err)
	}
	p := &loaded{
		listed: &listedPackage{ImportPath: "example.com/p"},
		fset:   fset, files: []*ast.File{file}, info: info, pkg: pkg,
	}
	var got []int
	for _, f := range cpuSpecificUses(p) {
		var line int
		if _, err := fmt.Sscanf(strings.TrimPrefix(f, "p.go:"), "%d", &line); err != nil {
			t.Fatalf("unparsable finding %q", f)
		}
		got = append(got, line)
	}
	sort.Ints(got)
	want := []int{10, 12, 14, 16, 17, 18}
	if fmt.Sprint(got) != fmt.Sprint(want) {
		t.Errorf("flagged lines %v, want %v", got, want)
	}
}
