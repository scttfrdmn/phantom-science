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

// TestNoFusedInstructions compiles internal/detmath and every output-path
// package for each architecture that fuses and fails if the object code
// contains a fused multiply-add. It
// checks the compiler's output, not only the source, so it also catches a
// fusion the source scan does not model. The standard library's math
// package is compiled the same way as a control, to show that the pattern
// matches the instructions the compiler really emits.
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
	for _, tg := range fusingTargets {
		name := tg.goarch + tg.goamd64
		t.Run(name, func(t *testing.T) {
			dump := func(pkg string) string {
				out := filepath.Join(tmp, name+"-"+strings.ReplaceAll(pkg, "/", "_")+".a")
				build := exec.Command(goBin, "build", "-o", out, pkg)
				build.Dir = root
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
			if n := len(fusedOp.FindAllString(dump("math"), -1)); n == 0 {
				t.Fatalf("control: no fused instructions found in package math for %s; the pattern %v no longer matches", name, fusedOp)
			}
			for _, pkg := range pkgs {
				for _, line := range strings.Split(dump(pkg), "\n") {
					if fusedOp.MatchString(line) {
						t.Errorf("fused instruction in %s: %s", pkg, strings.Join(strings.Fields(line), " "))
					}
				}
			}
		})
	}
}
