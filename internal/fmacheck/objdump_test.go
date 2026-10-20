//go:build !release

package fmacheck

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
)

// fusingTargets are the GOARCH settings on which the gc compiler emits fused
// multiply-add instructions for x*y + z (amd64 only from GOAMD64=v3).
var fusingTargets = []struct{ goarch, goamd64 string }{
	{"arm64", ""},
	{"amd64", "v3"},
	{"ppc64le", ""},
	{"s390x", ""},
	{"riscv64", ""},
	{"loong64", ""},
}

// fusedOp matches the fused multiply-add mnemonics of those architectures in
// `go tool objdump` output: FMADDD/FMSUBD/FNMADDD/FNMSUBD (arm64, riscv64,
// loong64), FMADD/FMSUB/FNMADD/FNMSUB (ppc64, s390x, and their single
// precision forms), and VFMADD231SD and friends (amd64).
var fusedOp = regexp.MustCompile(`\t(V?FN?M(ADD|SUB)[A-Z0-9]*)\s`)

// vexFMA matches the encoding of the amd64 VEX FMA instructions
// (VFMADD/VFMSUB/VFNMADD/VFNMSUB 132/213/231 and their packed forms: VEX
// prefix c4, map 0F38 (e2 with R=X=B=1, or any c4 xx with map 2), opcodes
// 96–9f, a6–af, b6–bf) in objdump's hex column. Some toolchains'
// disassemblers do not decode these instructions (go1.26.9 prints them as
// MOVL), so the bytes are checked as well as the mnemonic.
var vexFMA = regexp.MustCompile(`\t\t\tc4[0-9a-f][2a]..(9[6-9a-f]|a[6-9a-f]|b[6-9a-f])`)

// controlSrc is a function the compiler fuses on every target in
// fusingTargets; it shows the scan can see fused instructions there.
const controlSrc = "package fc\n\nfunc F(a, b, c float64) float64 { return a*b + c }\n"

// TestNoFusedInstructions compiles internal/detmath and every output-path
// package for each architecture that fuses and fails if the object code
// contains a fused multiply-add. It
// checks the compiler's output, not only the source, so it also catches a
// fusion the source scan does not model. A one-line a*b + c function,
// compiled for each target, is the control: it shows the patterns match the
// instructions the compiler really emits.
func TestNoFusedInstructions(t *testing.T) {
	if testing.Short() {
		t.Skip("cross-compiles for six architectures")
	}
	goBin, err := exec.LookPath("go")
	if err != nil {
		t.Skipf("go command not available: %v", err)
	}
	root := moduleRoot(t)
	tmp := t.TempDir()
	list := exec.Command(goBin, append([]string{"list"}, outputPackages...)...)
	list.Dir = root
	listed, err := list.Output()
	if err != nil {
		t.Fatalf("go list: %v", err)
	}
	pkgs := strings.Fields(string(listed))
	control := filepath.Join(tmp, "control")
	if err := os.MkdirAll(control, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(control, "go.mod"), []byte("module fc\n\ngo 1.24\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(control, "fc.go"), []byte(controlSrc), 0o644); err != nil {
		t.Fatal(err)
	}
	fused := func(line string) bool { return fusedOp.MatchString(line) || vexFMA.MatchString(line) }
	for _, tg := range fusingTargets {
		name := tg.goarch + tg.goamd64
		t.Run(name, func(t *testing.T) {
			dump := func(dir, pkg string) string {
				out := filepath.Join(tmp, name+"-"+strings.ReplaceAll(pkg, "/", "_")+".a")
				build := exec.Command(goBin, "build", "-o", out, pkg)
				build.Dir = dir
				build.Env = append(os.Environ(), "GOOS=linux", "GOARCH="+tg.goarch, "CGO_ENABLED=0")
				if tg.goamd64 != "" {
					build.Env = append(build.Env, "GOAMD64="+tg.goamd64)
				}
				if b, err := build.CombinedOutput(); err != nil {
					t.Fatalf("go build %s: %v\n%s", pkg, err, b)
				}
				var stdout, stderr bytes.Buffer
				od := exec.Command(goBin, "tool", "objdump", out)
				od.Stdout, od.Stderr = &stdout, &stderr
				if err := od.Run(); err != nil {
					t.Fatalf("go tool objdump %s: %v\n%s", pkg, err, stderr.String())
				}
				return stdout.String()
			}
			n := 0
			for _, line := range strings.Split(dump(control, "."), "\n") {
				if fused(line) {
					n++
				}
			}
			if n == 0 {
				t.Fatalf("control: a*b + c compiled for %s shows no fused instruction; the patterns %v / %v no longer match this toolchain's output", name, fusedOp, vexFMA)
			}
			for _, pkg := range pkgs {
				for _, line := range strings.Split(dump(root, pkg), "\n") {
					if fused(line) {
						t.Errorf("fused instruction in %s: %s", pkg, strings.Join(strings.Fields(line), " "))
					}
				}
			}
		})
	}
}
