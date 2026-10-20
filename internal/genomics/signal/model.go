package signal

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// MaxK is the largest supported k-mer size (4^12 = 16.7M levels). ONT's
// current models are 9-mers.
const MaxK = 12

// PoreModel is a k-mer pore model: the expected current level, in pA, of
// every k-mer, indexed by its 2-bit encoding (see KmerIndex).
//
// A PoreModel is read-only after loading and safe for concurrent use.
type PoreModel struct {
	K    int       // k-mer length
	Mean []float32 // expected level (pA), len 4^K
	Stdv []float32 // per-sample level spread (pA), len 4^K
	// Name is a free-form label (the file path when loaded from a file).
	Name string
	// Missing is the number of k-mers absent from the source table and
	// filled in under WithAllowMissing (0 for a complete table).
	Missing int
}

// NumKmers returns 4^K.
func (m *PoreModel) NumKmers() int { return 1 << (2 * m.K) }

// baseCode maps A/C/G/T (and U, either case) to 0..3 and everything else
// to -1.
var baseCode = func() (t [256]int8) {
	for i := range t {
		t[i] = -1
	}
	for _, c := range []struct {
		b byte
		v int8
	}{{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'U', 3}} {
		t[c.b] = c.v
		t[c.b+'a'-'A'] = c.v
	}
	return t
}()

// KmerIndex returns the 2-bit index of kmer (A=0, C=1, G=2, T/U=3, first
// base most significant, so indices follow lexicographic order AAA..., AAC,
// ...). It returns -1 if kmer contains a non-ACGTU byte.
func KmerIndex(kmer []byte) int {
	idx := 0
	for _, b := range kmer {
		c := baseCode[b]
		if c < 0 {
			return -1
		}
		idx = idx<<2 | int(c)
	}
	return idx
}

// KmerString is the inverse of KmerIndex for a k-mer of length k.
func KmerString(idx, k int) string {
	b := make([]byte, k)
	for i := k - 1; i >= 0; i-- {
		b[i] = "ACGT"[idx&3]
		idx >>= 2
	}
	return string(b)
}

// LoadOption configures LoadPoreModel / ReadPoreModel.
type LoadOption func(*loadConfig)

type loadConfig struct {
	scale, shift float64
	defaultStdv  float64
	allowMissing bool
	expectK      int
}

// WithScaling converts levels from the table's units to pA:
// pA = level*scale + shift, stdv_pA = stdv*scale. Use it for ONT's
// kmer_models tables, which are in standardised units (about zero mean,
// unit spread); Profile.ModelOptions supplies preset values. Without it the
// table is assumed to be in pA already (squigulator's format).
func WithScaling(scale, shift float64) LoadOption {
	return func(c *loadConfig) { c.scale, c.shift = scale, shift }
}

// WithDefaultStdv sets the level stdv, in the table's units (before
// scaling), used for rows without a level_stdv column. The default is
// DefaultModelStdv.
func WithDefaultStdv(sd float64) LoadOption {
	return func(c *loadConfig) { c.defaultStdv = sd }
}

// WithAllowMissing accepts a table with fewer than 4^k k-mers; missing
// k-mers get the mean of the present levels and the median present stdv,
// and PoreModel.Missing records how many were filled. Without it a short
// table is an error.
func WithAllowMissing() LoadOption {
	return func(c *loadConfig) { c.allowMissing = true }
}

// WithK requires the table's k-mer length to be k.
func WithK(k int) LoadOption {
	return func(c *loadConfig) { c.expectK = k }
}

// DefaultModelStdv is the level stdv assigned to rows that have none, in
// table units. ONT's kmer_models tables carry only level means; 0.08
// standardised units is about 1.9 pA with the R10.4.1 scaling, in line with
// the 1.5–2.5 pA within-event spread of R10.4.1 signal.
const DefaultModelStdv = 0.08

// LoadPoreModel reads a k-mer level table from path (gzip-compressed if the
// name ends in .gz). See ReadPoreModel for the accepted formats.
func LoadPoreModel(path string, opts ...LoadOption) (*PoreModel, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied model path
	if err != nil {
		return nil, fmt.Errorf("pore model: %w", err)
	}
	defer f.Close()
	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("pore model %s: %w", path, err)
		}
		defer gz.Close()
		r = gz
	}
	m, err := ReadPoreModel(r, opts...)
	if err != nil {
		return nil, fmt.Errorf("pore model %s: %w", path, err)
	}
	m.Name = path
	return m, nil
}

// ReadPoreModel parses a k-mer level table. Accepted formats:
//
//   - ONT kmer_models (github.com/nanoporetech/kmer_models): no header,
//     two whitespace-separated columns "kmer level" in standardised units;
//   - squigulator / nanopolish style: optional "#..." metadata lines
//     (a "#k N" line is checked against the k-mer length), an optional
//     header row starting with "kmer", and columns
//     "kmer level_mean level_stdv [extra columns ignored]", in pA.
//
// k is inferred from the first k-mer; every k-mer must have that length,
// consist of A/C/G/T (U is read as T), and appear once. All 4^k k-mers are
// required unless WithAllowMissing is given. Rows may be in any order.
func ReadPoreModel(r io.Reader, opts ...LoadOption) (*PoreModel, error) {
	cfg := loadConfig{scale: 1, defaultStdv: DefaultModelStdv}
	for _, o := range opts {
		o(&cfg)
	}
	if !(cfg.scale > 0) || math.IsInf(cfg.scale, 0) || math.IsNaN(cfg.shift) || math.IsInf(cfg.shift, 0) {
		return nil, fmt.Errorf("invalid scaling scale=%g shift=%g", cfg.scale, cfg.shift)
	}
	if !(cfg.defaultStdv >= 0) {
		return nil, fmt.Errorf("invalid default stdv %g", cfg.defaultStdv)
	}

	var (
		m        *PoreModel
		seen     []bool
		n        int
		headerK  int
		lineNo   int
		dataSeen bool
	)
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 1024*1024)
	for sc.Scan() {
		lineNo++
		line := strings.TrimSpace(sc.Text())
		if line == "" {
			continue
		}
		if line[0] == '#' {
			f := strings.Fields(line)
			if len(f) >= 2 && f[0] == "#k" {
				k, err := strconv.Atoi(f[1])
				if err != nil || k <= 0 {
					return nil, fmt.Errorf("line %d: invalid #k %q", lineNo, f[1])
				}
				headerK = k
			}
			continue
		}
		f := strings.Fields(line)
		if !dataSeen && strings.EqualFold(f[0], "kmer") {
			continue // header row
		}
		dataSeen = true
		if len(f) < 2 {
			return nil, fmt.Errorf("line %d: want at least 2 columns (kmer level), got %d", lineNo, len(f))
		}
		if m == nil {
			k := len(f[0])
			if k < 1 || k > MaxK {
				return nil, fmt.Errorf("line %d: k-mer length %d outside 1..%d", lineNo, k, MaxK)
			}
			if headerK != 0 && headerK != k {
				return nil, fmt.Errorf("line %d: k-mer length %d disagrees with #k %d", lineNo, k, headerK)
			}
			if cfg.expectK != 0 && cfg.expectK != k {
				return nil, fmt.Errorf("line %d: k-mer length %d, want k=%d", lineNo, k, cfg.expectK)
			}
			nk := 1 << (2 * k)
			m = &PoreModel{K: k, Mean: make([]float32, nk), Stdv: make([]float32, nk)}
			seen = make([]bool, nk)
		}
		if len(f[0]) != m.K {
			return nil, fmt.Errorf("line %d: k-mer %q has length %d, want %d", lineNo, f[0], len(f[0]), m.K)
		}
		idx := KmerIndex([]byte(f[0]))
		if idx < 0 {
			return nil, fmt.Errorf("line %d: k-mer %q contains a non-ACGT base", lineNo, f[0])
		}
		if seen[idx] {
			return nil, fmt.Errorf("line %d: duplicate k-mer %s", lineNo, f[0])
		}
		mean, err := parseFinite(f[1])
		if err != nil {
			return nil, fmt.Errorf("line %d: level_mean: %w", lineNo, err)
		}
		sd := cfg.defaultStdv
		if len(f) >= 3 {
			if sd, err = parseFinite(f[2]); err != nil {
				return nil, fmt.Errorf("line %d: level_stdv: %w", lineNo, err)
			}
			if sd < 0 {
				return nil, fmt.Errorf("line %d: negative level_stdv %g", lineNo, sd)
			}
		}
		seen[idx] = true
		m.Mean[idx] = float32(float64(mean*cfg.scale) + cfg.shift)
		m.Stdv[idx] = float32(sd * cfg.scale)
		n++
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if m == nil {
		return nil, fmt.Errorf("no k-mer rows")
	}
	if missing := m.NumKmers() - n; missing > 0 {
		if !cfg.allowMissing {
			return nil, fmt.Errorf("%d of %d %d-mers present (%d missing); a complete table is required", n, m.NumKmers(), m.K, missing)
		}
		var sum float64
		sds := make([]float64, 0, n)
		for i, ok := range seen {
			if ok {
				sum += float64(m.Mean[i])
				sds = append(sds, float64(m.Stdv[i]))
			}
		}
		sort.Float64s(sds)
		fillMean, fillSD := float32(sum/float64(n)), float32(sds[len(sds)/2])
		for i, ok := range seen {
			if !ok {
				m.Mean[i], m.Stdv[i] = fillMean, fillSD
			}
		}
		m.Missing = missing
	}
	return m, nil
}

func parseFinite(s string) (float64, error) {
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, fmt.Errorf("%q is not a number", s)
	}
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return 0, fmt.Errorf("%q is not finite", s)
	}
	return v, nil
}
