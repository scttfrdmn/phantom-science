//go:build !release

package fmacheck

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"go/ast"
	"go/format"
	"go/importer"
	"go/parser"
	"go/token"
	"go/types"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// outputPackages are the packages whose floating-point arithmetic affects
// generated output (sequences, qualities, read lengths, signal, expression,
// genotypes). Keep this list in step with the generation path.
var outputPackages = []string{
	"./internal/genomics",
	"./internal/genomics/authenticity",
	"./internal/genomics/clinical",
	"./internal/genomics/illumina",
	"./internal/genomics/nanopore/...",
	"./internal/genomics/population",
	"./internal/genomics/readlen",
	"./internal/genomics/reference",
	"./internal/genomics/seqread",
	"./internal/genomics/signal",
	"./internal/genomics/slow5",
	"./internal/genomics/transcriptome",
	"./internal/population",
}

var fix = flag.Bool("fix", false, "rewrite flagged products as explicit conversions in place")

// TestNoFusableFloatProducts fails for every floating-point product on the
// output path that the compiler could fuse into a multiply-add, which would
// make output depend on the CPU architecture. Wrap the product in an
// explicit conversion (float64(a*b) + c) to forbid the fusion.
//
// Run with -args -fix to apply the conversions automatically; review and
// gofmt the result.
func TestNoFusableFloatProducts(t *testing.T) {
	pkgs := loadPackages(t, outputPackages)
	var all []finding
	for _, p := range pkgs {
		all = append(all, checkPackage(t, p)...)
	}
	if *fix {
		applyFixes(t, all)
		return
	}
	for _, f := range all {
		t.Errorf("%s: %s: %s; wrap the product in %s(...)", f.pos, f.why, f.text, f.conv)
	}
	if len(all) > 0 {
		t.Logf("%d fusable floating-point products; see internal/fmacheck/doc.go", len(all))
	}
}

// listedPackage is the subset of `go list -json` output used here.
type listedPackage struct {
	ImportPath string
	Dir        string
	Export     string
	GoFiles    []string
	CgoFiles   []string
	DepOnly    bool
	Module     *struct{ Path string }
	Error      *struct{ Err string }
}

type loaded struct {
	listed *listedPackage
	fset   *token.FileSet
	files  []*ast.File
	info   *types.Info
	pkg    *types.Package
}

func moduleRoot(t *testing.T) string {
	t.Helper()
	dir, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatal("go.mod not found")
		}
		dir = parent
	}
}

func loadPackages(t *testing.T, patterns []string) []*loaded {
	t.Helper()
	goBin, err := exec.LookPath("go")
	if err != nil {
		t.Skipf("go command not available: %v", err)
	}
	root := moduleRoot(t)
	args := append([]string{"list", "-e", "-export", "-deps", "-json"}, patterns...)
	cmd := exec.Command(goBin, args...)
	cmd.Dir = root
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("go list: %v\n%s", err, stderr.String())
	}

	exports := map[string]string{}
	var targets []*listedPackage
	dec := json.NewDecoder(bytes.NewReader(out))
	for {
		var p listedPackage
		if err := dec.Decode(&p); err == io.EOF {
			break
		} else if err != nil {
			t.Fatalf("decoding go list output: %v", err)
		}
		if p.Error != nil {
			t.Fatalf("go list %s: %s", p.ImportPath, p.Error.Err)
		}
		if p.Export != "" {
			exports[p.ImportPath] = p.Export
		}
		if !p.DepOnly {
			pp := p
			targets = append(targets, &pp)
		}
	}
	if len(targets) == 0 {
		t.Fatal("no packages matched")
	}

	fset := token.NewFileSet()
	imp := importer.ForCompiler(fset, "gc", func(path string) (io.ReadCloser, error) {
		file, ok := exports[path]
		if !ok {
			return nil, fmt.Errorf("no export data for %q", path)
		}
		return os.Open(file)
	})

	var res []*loaded
	for _, lp := range targets {
		var files []*ast.File
		for _, name := range append(append([]string{}, lp.GoFiles...), lp.CgoFiles...) {
			f, err := parser.ParseFile(fset, filepath.Join(lp.Dir, name), nil, parser.ParseComments)
			if err != nil {
				t.Fatal(err)
			}
			files = append(files, f)
		}
		info := &types.Info{
			Types: map[ast.Expr]types.TypeAndValue{},
			Defs:  map[*ast.Ident]types.Object{},
			Uses:  map[*ast.Ident]types.Object{},
		}
		conf := types.Config{Importer: imp, FakeImportC: true}
		pkg, err := conf.Check(lp.ImportPath, fset, files, info)
		if err != nil {
			t.Fatalf("type-checking %s: %v", lp.ImportPath, err)
		}
		res = append(res, &loaded{listed: lp, fset: fset, files: files, info: info, pkg: pkg})
	}
	return res
}

type finding struct {
	pos  token.Position
	why  string
	text string
	conv string
	// Byte offsets of the expression to wrap and, when the product is
	// directly parenthesised, the offset of the outermost '('.
	file       string
	start, end int
	lparen     int
}

// isFloat reports whether t is a (possibly named) float32 or float64.
func isFloat(t types.Type) bool {
	b, ok := t.Underlying().(*types.Basic)
	return ok && b.Info()&types.IsFloat != 0 && b.Info()&types.IsUntyped == 0
}

func checkPackage(t *testing.T, p *loaded) []finding {
	t.Helper()
	modPath := ""
	if p.listed.Module != nil {
		modPath = p.listed.Module.Path
	}
	var out []finding
	for _, file := range p.files {
		c := &checker{p: p, modPath: modPath, file: file}
		c.collectUses()
		var stack []ast.Node
		ast.Inspect(file, func(n ast.Node) bool {
			if n == nil {
				stack = stack[:len(stack)-1]
				return false
			}
			stack = append(stack, n)
			c.visit(n, stack)
			return true
		})
		out = append(out, c.out...)
	}
	return out
}

type checker struct {
	p       *loaded
	modPath string
	file    *ast.File
	// uses maps a local variable to the identifiers that read it, together
	// with each identifier's ancestor chain.
	uses map[*types.Var][][]ast.Node
	out  []finding
}

func (c *checker) collectUses() {
	c.uses = map[*types.Var][][]ast.Node{}
	var stack []ast.Node
	ast.Inspect(c.file, func(n ast.Node) bool {
		if n == nil {
			stack = stack[:len(stack)-1]
			return false
		}
		stack = append(stack, n)
		if id, ok := n.(*ast.Ident); ok {
			if v, ok := c.p.info.Uses[id].(*types.Var); ok && !v.IsField() {
				c.uses[v] = append(c.uses[v], append([]ast.Node(nil), stack...))
			}
		}
		return true
	})
}

func (c *checker) isConst(e ast.Expr) bool {
	tv, ok := c.p.info.Types[e]
	return ok && tv.Value != nil
}

func (c *checker) floatType(e ast.Expr) (types.Type, bool) {
	tv, ok := c.p.info.Types[e]
	if !ok || tv.Type == nil || !isFloat(tv.Type) {
		return nil, false
	}
	return tv.Type, true
}

func (c *checker) visit(n ast.Node, stack []ast.Node) {
	switch e := n.(type) {
	case *ast.BinaryExpr:
		if e.Op != token.MUL || c.isConst(e) {
			return
		}
		typ, ok := c.floatType(e)
		if !ok {
			return
		}
		if why := c.unsafeConsumer(stack); why != "" {
			c.report(e, typ, why, stack)
		}
	case *ast.AssignStmt:
		// x *= y stores an unrounded product in x.
		if e.Tok != token.MUL_ASSIGN || len(e.Lhs) != 1 {
			return
		}
		typ, ok := c.floatType(e.Lhs[0])
		if !ok {
			return
		}
		if c.storeIsUnsafe(e.Lhs[0]) {
			pos := c.p.fset.Position(e.Pos())
			c.out = append(c.out, finding{
				pos:  pos,
				why:  "product stored by *= may later be fused",
				text: c.src(e),
				conv: c.convName(typ),
			})
		}
	}
}

// strip walks up from the product past parentheses and unary +/-, which do
// not stop fusion, and returns the consuming node and the child it consumes.
func strip(stack []ast.Node) (parent, child ast.Node, i int) {
	i = len(stack) - 1
	child = stack[i]
	for i > 0 {
		p := stack[i-1]
		switch pe := p.(type) {
		case *ast.ParenExpr:
			child, i = p, i-1
			continue
		case *ast.UnaryExpr:
			if pe.Op == token.SUB || pe.Op == token.ADD {
				child, i = p, i-1
				continue
			}
		}
		return p, child, i - 1
	}
	return nil, child, -1
}

// unsafeConsumer returns why the product at the top of stack could be fused
// into an add or subtract, or "" if it cannot.
func (c *checker) unsafeConsumer(stack []ast.Node) string {
	parent, child, _ := strip(stack)
	switch p := parent.(type) {
	case *ast.BinaryExpr:
		switch p.Op {
		case token.ADD, token.SUB:
			return "product is an operand of " + p.Op.String()
		}
		// Operand of *, /, or a comparison: the product is consumed by an
		// operation that cannot absorb it.
		return ""
	case *ast.CallExpr:
		if c.isConversion(p) {
			return ""
		}
		if c.isSafeCallee(p) {
			return ""
		}
		return "product is an argument to a call that may be inlined"
	case *ast.AssignStmt:
		switch p.Tok {
		case token.ADD_ASSIGN, token.SUB_ASSIGN:
			for _, r := range p.Rhs {
				if r == child {
					return "product is the operand of " + p.Tok.String()
				}
			}
			return ""
		case token.MUL_ASSIGN, token.QUO_ASSIGN:
			return ""
		case token.ASSIGN, token.DEFINE:
			for k, r := range p.Rhs {
				if r == child && len(p.Lhs) == len(p.Rhs) {
					if c.storeIsUnsafe(p.Lhs[k]) {
						return "product is stored and later added or subtracted"
					}
					return ""
				}
			}
		}
		return "product is stored"
	case *ast.ValueSpec:
		for k, v := range p.Values {
			if v == child && k < len(p.Names) {
				if c.storeIsUnsafe(p.Names[k]) {
					return "product is stored and later added or subtracted"
				}
				return ""
			}
		}
		return "product is stored"
	case *ast.IndexExpr, *ast.IndexListExpr, *ast.SliceExpr, *ast.ExprStmt:
		return ""
	case *ast.SwitchStmt, *ast.IfStmt, *ast.ForStmt:
		return ""
	case nil:
		return ""
	}
	// Returns, composite literals, sends, field stores: the value escapes
	// the expression and may reach an add after inlining or store-to-load
	// forwarding.
	return "product escapes the expression"
}

func (c *checker) isConversion(call *ast.CallExpr) bool {
	tv, ok := c.p.info.Types[call.Fun]
	return ok && tv.IsType()
}

// isSafeCallee reports whether a product passed to call cannot be fused by
// inlining: functions outside this module (the standard library's float
// helpers are intrinsics or too large to inline into an add of their
// argument) and fmt-style variadic interfaces.
func (c *checker) isSafeCallee(call *ast.CallExpr) bool {
	var obj types.Object
	switch f := ast.Unparen(call.Fun).(type) {
	case *ast.Ident:
		obj = c.p.info.Uses[f]
	case *ast.SelectorExpr:
		obj = c.p.info.Uses[f.Sel]
	}
	fn, ok := obj.(*types.Func)
	if !ok {
		if b, ok := obj.(*types.Builtin); ok {
			// min/max are exact; append/copy store values.
			return b.Name() == "min" || b.Name() == "max"
		}
		return false
	}
	if fn.Pkg() == nil {
		return true
	}
	path := fn.Pkg().Path()
	return c.modPath != "" && path != c.modPath && !strings.HasPrefix(path, c.modPath+"/")
}

// storeIsUnsafe reports whether storing an unrounded product into lhs could
// let the compiler fuse it into a later add or subtract. Only local
// variables whose every read is a safe consumer are considered safe.
func (c *checker) storeIsUnsafe(lhs ast.Expr) bool {
	id, ok := ast.Unparen(lhs).(*ast.Ident)
	if !ok {
		return true
	}
	if id.Name == "_" {
		return false
	}
	obj := c.p.info.Defs[id]
	if obj == nil {
		obj = c.p.info.Uses[id]
	}
	v, ok := obj.(*types.Var)
	if !ok || v.IsField() || v.Parent() == nil || v.Parent() == c.p.pkg.Scope() {
		return true
	}
	for _, st := range c.uses[v] {
		if c.readIsUnsafe(st) {
			return true
		}
	}
	return false
}

// readIsUnsafe reports whether the read of a variable at the top of stack
// could absorb a product held in that variable.
func (c *checker) readIsUnsafe(stack []ast.Node) bool {
	parent, child, _ := strip(stack)
	switch p := parent.(type) {
	case *ast.BinaryExpr:
		return p.Op == token.ADD || p.Op == token.SUB
	case *ast.CallExpr:
		if p.Fun == child {
			return false
		}
		return !c.isConversion(p) && !c.isSafeCallee(p)
	case *ast.AssignStmt:
		switch p.Tok {
		case token.ADD_ASSIGN, token.SUB_ASSIGN:
			return true // both x += c with x a product and c += x
		case token.MUL_ASSIGN, token.QUO_ASSIGN:
			return false
		}
		for _, l := range p.Lhs {
			if l == child {
				return false // overwritten, not read
			}
		}
		return true // copied elsewhere
	case *ast.IncDecStmt:
		return true
	case *ast.IndexExpr, *ast.IndexListExpr, *ast.SliceExpr, *ast.ExprStmt,
		*ast.SwitchStmt, *ast.IfStmt, *ast.ForStmt, *ast.SelectorExpr:
		return false
	}
	return true
}

func (c *checker) convName(t types.Type) string {
	return types.TypeString(t, types.RelativeTo(c.p.pkg))
}

func (c *checker) src(n ast.Node) string {
	var b bytes.Buffer
	if err := format.Node(&b, c.p.fset, n); err != nil {
		return "?"
	}
	s := b.String()
	if len(s) > 80 {
		s = s[:77] + "..."
	}
	return s
}

func (c *checker) report(e *ast.BinaryExpr, typ types.Type, why string, stack []ast.Node) {
	pos := c.p.fset.Position(e.Pos())
	f := finding{
		pos:    pos,
		why:    why,
		text:   c.src(e),
		conv:   c.convName(typ),
		file:   pos.Filename,
		start:  pos.Offset,
		end:    c.p.fset.Position(e.End()).Offset,
		lparen: -1,
	}
	// If the product is directly parenthesised, turn the outermost
	// parentheses into the conversion instead of nesting them.
	for i := len(stack) - 2; i >= 0; i-- {
		pe, ok := stack[i].(*ast.ParenExpr)
		if !ok {
			break
		}
		f.lparen = c.p.fset.Position(pe.Lparen).Offset
	}
	c.out = append(c.out, f)
}

type edit struct {
	off  int
	text string
}

func applyFixes(t *testing.T, all []finding) {
	t.Helper()
	byFile := map[string][]edit{}
	for _, f := range all {
		if f.file == "" {
			t.Logf("%s: not fixed automatically: %s", f.pos, f.text)
			continue
		}
		if f.lparen >= 0 {
			byFile[f.file] = append(byFile[f.file], edit{f.lparen, f.conv})
			continue
		}
		byFile[f.file] = append(byFile[f.file],
			edit{f.start, f.conv + "("}, edit{f.end, ")"})
	}
	for file, edits := range byFile {
		src, err := os.ReadFile(file)
		if err != nil {
			t.Fatal(err)
		}
		sort.SliceStable(edits, func(i, j int) bool { return edits[i].off > edits[j].off })
		for _, e := range edits {
			src = append(src[:e.off], append([]byte(e.text), src[e.off:]...)...)
		}
		if out, err := format.Source(src); err == nil {
			src = out
		} else {
			t.Errorf("%s: gofmt after fix: %v", file, err)
		}
		if err := os.WriteFile(file, src, 0o644); err != nil {
			t.Fatal(err)
		}
		t.Logf("fixed %d sites in %s", len(edits), file)
	}
}

// TestCheckerFindsFusableSites runs the checker on a small source file with
// known fusable and safe products.
func TestCheckerFindsFusableSites(t *testing.T) {
	const src = `package p

import "math"

type pA float64

func helper(x float64) float64 { return x }

func f(a, b, c float64, s []float64, q pA) (float64, float64) {
	r1 := a*b + c          // flag
	r2 := c - a*b*c        // flag
	r3 := float64(a*b) + c // safe
	r4 := (a * b) / c      // safe
	r5 := -(a * b) + c     // flag
	var acc float64
	acc += a * b          // flag
	acc -= float64(a * b) // safe
	t := a * b            // flag: t is added below
	r6 := t + c
	u := a * b // safe: u is only compared
	ok := u < c
	s[0] *= a              // flag: stored element
	r7 := math.Sqrt(a * b) // safe: math call
	r8 := helper(a * b)    // flag: may be inlined
	q2 := q*q + q          // flag: named float type
	const k = 2.0 * 3.0
	r9 := k + c   // safe: constant product
	i := 3*4 + 1  // safe: integer
	_ = ok
	_, _, _, _, _, _, _, _, _, _, _ = r1, r2, r3, r4, r5, r7, r8, q2, r9, i, acc
	return r6, a * b // flag: returned
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
		listed: &listedPackage{ImportPath: "example.com/p", Module: &struct{ Path string }{"example.com/p"}},
		fset:   fset, files: []*ast.File{file}, info: info, pkg: pkg,
	}
	var got []int
	for _, f := range checkPackage(t, p) {
		got = append(got, f.pos.Line)
	}
	want := []int{10, 11, 14, 16, 18, 22, 24, 25, 31}
	if fmt.Sprint(got) != fmt.Sprint(want) {
		t.Errorf("flagged lines %v, want %v", got, want)
	}
}
