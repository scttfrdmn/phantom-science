//go:build !release

package signal

import (
	"compress/gzip"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// testLevel is the deterministic synthetic level table used by the tests:
// standardised-looking values in roughly [-2, 2].
func testLevel(idx int) float64 {
	return 4*(math.Mod(float64(idx)*0.6180339887, 1)) - 2
}

// testModelText renders a complete k-mer table. cols is 2 (ONT) or 3
// (with stdv); header adds squigulator-style metadata and a header row.
func testModelText(k, cols int, header bool) string {
	var b strings.Builder
	if header {
		fmt.Fprintf(&b, "#ont_model_name\ttest\n#k\t%d\nkmer\tlevel_mean\tlevel_stdv\tsd_mean\tsd_stdv\n", k)
	}
	for i := 0; i < 1<<(2*k); i++ {
		km := KmerString(i, k)
		switch cols {
		case 2:
			fmt.Fprintf(&b, "%s\t%.6f\n", km, testLevel(i))
		default:
			fmt.Fprintf(&b, "%s\t%.6f\t%.6f\t1.0\t0.1\n", km, 85+20*testLevel(i), 1+0.3*float64(i%5))
		}
	}
	return b.String()
}

func writeFile(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o600); err != nil {
		t.Fatal(err)
	}
	return p
}

// testModel is a k=5 model in pA: levels 45–125 pA, stdv 1.0–2.2 pA.
func testModel(t testing.TB) *PoreModel {
	t.Helper()
	m, err := ReadPoreModel(strings.NewReader(testModelText(5, 3, true)))
	if err != nil {
		t.Fatal(err)
	}
	return m
}

func TestKmerIndexRoundTrip(t *testing.T) {
	for k := 1; k <= 6; k++ {
		for i := 0; i < 1<<(2*k); i += 7 {
			if got := KmerIndex([]byte(KmerString(i, k))); got != i {
				t.Fatalf("k=%d: KmerIndex(KmerString(%d)) = %d", k, i, got)
			}
		}
	}
	if KmerIndex([]byte("acgu")) != KmerIndex([]byte("ACGT")) {
		t.Error("lowercase / U not equivalent to ACGT")
	}
	if KmerIndex([]byte("ACNT")) != -1 {
		t.Error("N should give -1")
	}
	// Lexicographic order: the ONT tables list AAAA.., AAAC.., so index 1 is ...AC.
	if KmerString(1, 3) != "AAC" || KmerString(63, 3) != "TTT" {
		t.Errorf("KmerString order wrong: %s %s", KmerString(1, 3), KmerString(63, 3))
	}
}

func TestLoadPoreModelONTFormatScaling(t *testing.T) {
	path := writeFile(t, "9mer_levels_v1.txt", testModelText(4, 2, false))
	const scale, shift = 23.5, 93.7
	m, err := LoadPoreModel(path, WithScaling(scale, shift))
	if err != nil {
		t.Fatal(err)
	}
	if m.K != 4 || len(m.Mean) != 256 || m.Missing != 0 {
		t.Fatalf("K=%d len=%d missing=%d", m.K, len(m.Mean), m.Missing)
	}
	for _, i := range []int{0, 1, 100, 255} {
		want := testLevel(i)*scale + shift
		if math.Abs(float64(m.Mean[i])-want) > 1e-3 {
			t.Errorf("mean[%d] = %g, want %g", i, m.Mean[i], want)
		}
		if math.Abs(float64(m.Stdv[i])-DefaultModelStdv*scale) > 1e-5 {
			t.Errorf("stdv[%d] = %g, want default %g", i, m.Stdv[i], DefaultModelStdv*scale)
		}
	}
	m2, err := LoadPoreModel(path, WithScaling(scale, shift), WithDefaultStdv(0.1))
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(float64(m2.Stdv[3])-0.1*scale) > 1e-5 {
		t.Errorf("WithDefaultStdv: stdv = %g, want %g", m2.Stdv[3], 0.1*scale)
	}
	// Presets convert ONT tables via ModelOptions (k checked).
	p, _ := Preset(PresetDNAR10)
	if _, err := LoadPoreModel(path, p.ModelOptions()...); err == nil || !strings.Contains(err.Error(), "want k=9") {
		t.Errorf("preset ModelOptions should reject a 4-mer table, got %v", err)
	}
}

func TestLoadPoreModelSquigulatorFormatAndGzip(t *testing.T) {
	text := testModelText(3, 3, true)
	m, err := ReadPoreModel(strings.NewReader(text))
	if err != nil {
		t.Fatal(err)
	}
	if m.K != 3 {
		t.Fatalf("K = %d", m.K)
	}
	idx := KmerIndex([]byte("CGT"))
	if want := 85 + 20*testLevel(idx); math.Abs(float64(m.Mean[idx])-want) > 1e-3 {
		t.Errorf("CGT mean %g want %g (pA, unscaled)", m.Mean[idx], want)
	}
	if want := 1 + 0.3*float64(idx%5); math.Abs(float64(m.Stdv[idx])-want) > 1e-5 {
		t.Errorf("CGT stdv %g want %g", m.Stdv[idx], want)
	}

	// Rows in any order, CRLF and space separation.
	lines := strings.Split(strings.TrimSpace(testModelText(3, 2, false)), "\n")
	for i, j := 0, len(lines)-1; i < j; i, j = i+1, j-1 {
		lines[i], lines[j] = lines[j], lines[i]
	}
	rev := strings.ReplaceAll(strings.Join(lines, "\r\n"), "\t", "  ")
	m2, err := ReadPoreModel(strings.NewReader(rev))
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(float64(m2.Mean[idx])-testLevel(idx)) > 1e-5 {
		t.Errorf("reversed order: mean %g want %g", m2.Mean[idx], testLevel(idx))
	}

	gzPath := filepath.Join(t.TempDir(), "model.txt.gz")
	f, err := os.Create(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	zw := gzip.NewWriter(f)
	_, _ = zw.Write([]byte(text))
	_ = zw.Close()
	_ = f.Close()
	mg, err := LoadPoreModel(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	if mg.Mean[idx] != m.Mean[idx] || mg.Name != gzPath {
		t.Errorf("gzip model differs: %g vs %g (name %q)", mg.Mean[idx], m.Mean[idx], mg.Name)
	}
}

func TestLoadPoreModelErrors(t *testing.T) {
	full := testModelText(2, 2, false)
	rows := strings.Split(strings.TrimSpace(full), "\n")
	cases := map[string]struct {
		text string
		opts []LoadOption
		want string
	}{
		"missing rows":     {strings.Join(rows[:15], "\n"), nil, "15 of 16"},
		"duplicate":        {full + rows[0] + "\n", nil, "duplicate"},
		"bad level":        {strings.Replace(full, rows[3], "ATG\tx", 1), nil, "length 3"},
		"not a number":     {strings.Replace(full, rows[3], "AT\tabc", 1), nil, "not a number"},
		"nan level":        {strings.Replace(full, rows[3], "AT\tNaN", 1), nil, "not finite"},
		"one column":       {strings.Replace(full, rows[3], "AT", 1), nil, "at least 2 columns"},
		"non ACGT":         {strings.Replace(full, rows[3], "AN\t1.0", 1), nil, "non-ACGT"},
		"negative stdv":    {strings.Replace(full, rows[3], "AT\t1.0\t-1", 1), nil, "negative"},
		"k mismatch":       {"#k\t3\n" + full, nil, "disagrees with #k"},
		"WithK mismatch":   {full, []LoadOption{WithK(9)}, "want k=9"},
		"empty":            {"# nothing\n\n", nil, "no k-mer rows"},
		"bad scaling":      {full, []LoadOption{WithScaling(0, 1)}, "invalid scaling"},
		"too long k":       {strings.Repeat("A", MaxK+1) + "\t1\n", nil, "outside"},
		"header after row": {rows[0] + "\nkmer\tlevel\n" + strings.Join(rows[1:], "\n"), nil, "length 4"},
	}
	for name, c := range cases {
		t.Run(name, func(t *testing.T) {
			_, err := ReadPoreModel(strings.NewReader(c.text), c.opts...)
			if err == nil || !strings.Contains(err.Error(), c.want) {
				t.Fatalf("err = %v, want containing %q", err, c.want)
			}
		})
	}
	if _, err := LoadPoreModel(filepath.Join(t.TempDir(), "absent.txt")); err == nil {
		t.Error("missing file: want error")
	}
}

func TestLoadPoreModelAllowMissing(t *testing.T) {
	rows := strings.Split(strings.TrimSpace(testModelText(2, 3, false)), "\n")
	m, err := ReadPoreModel(strings.NewReader(strings.Join(rows[:12], "\n")), WithAllowMissing())
	if err != nil {
		t.Fatal(err)
	}
	if m.Missing != 4 {
		t.Fatalf("Missing = %d, want 4", m.Missing)
	}
	var sum float64
	for i := 0; i < 12; i++ {
		sum += float64(m.Mean[i])
	}
	for i := 12; i < 16; i++ {
		if math.Abs(float64(m.Mean[i])-sum/12) > 1e-3 || m.Stdv[i] <= 0 {
			t.Errorf("filled k-mer %d: mean %g (want %g) stdv %g", i, m.Mean[i], sum/12, m.Stdv[i])
		}
	}
}
