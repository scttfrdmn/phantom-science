//go:build !release && cgo
// +build !release,cgo

package fingerprinting

import (
	"errors"
	"fmt"
	"runtime/debug"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

const steganoModule = "github.com/scttfrdmn/stegano/go"

var _ SAMFingerprintingEngine = (*SteganoEngine)(nil)

// fastqRecords builds n FASTQ records of length l whose Phred+33 qualities
// come from qual(read, position).
func fastqRecords(n, l int, qual func(read, pos int) byte) []byte {
	var b strings.Builder
	for i := 0; i < n; i++ {
		seq := make([]byte, l)
		q := make([]byte, l)
		for j := 0; j < l; j++ {
			seq[j] = "ACGT"[(i*7+j*3)%4]
			q[j] = qual(i, j)
		}
		fmt.Fprintf(&b, "@read_%d 1:N:0:ACGT\n%s\n+\n%s\n", i, seq, q)
	}
	return []byte(b.String())
}

func newTestEngine(t *testing.T) *SteganoEngine {
	t.Helper()
	e, err := NewSteganoEngine("customer-123")
	require.NoError(t, err)
	t.Cleanup(e.Close)
	return e
}

func requireCode(t *testing.T, err error, code ErrorCode) *FingerprintingError {
	t.Helper()
	var fe *FingerprintingError
	require.True(t, errors.As(err, &fe), "want *FingerprintingError, got %T: %v", err, err)
	assert.Equal(t, code, fe.Code, "error: %v", err)
	return fe
}

// The Go module and the libstegano it links against must come from the same
// release, or the C header and library symbols disagree. A stale library left
// in a shared cargo target directory is the usual cause.
func TestLinkedLibraryMatchesGoModule(t *testing.T) {
	info, ok := debug.ReadBuildInfo()
	require.True(t, ok, "build info unavailable")
	var modVersion string
	for _, dep := range info.Deps {
		if dep.Path == steganoModule {
			modVersion = dep.Version
		}
	}
	require.NotEmpty(t, modVersion, "%s not found in build info", steganoModule)
	assert.Equal(t, modVersion, "v"+LibraryVersion(),
		"linked libstegano does not match the Go module version; rebuild it with scripts/build-stegano.sh")
}

// With CGO, New must return the real engine. Falling back to the legacy engine
// would silently produce output with no extractable fingerprint.
func TestNew_UsesSteganoEngine(t *testing.T) {
	e := New("customer-123")
	defer e.Close()
	assert.IsType(t, &SteganoEngine{}, e)
}

func TestSteganoEngine_GenerateFingerprint(t *testing.T) {
	fp, err := newTestEngine(t).GenerateFingerprint("dataset-001")
	require.NoError(t, err)
	assert.True(t, strings.HasPrefix(fp, "customer-123:dataset-001:stegano-"+LibraryVersion()+":ts-"), fp)
}

func TestSteganoEngine_FASTQRoundTrip(t *testing.T) {
	cases := []struct {
		name string
		data []byte
	}{
		// NovaSeq-like: 150 bp, high quality.
		{"novaseq_150bp", fastqRecords(50, 150, func(_, j int) byte { return byte('F' - j%3) })},
		// MiSeq-like: 300 bp, quality declining along the read.
		{"miseq_300bp", fastqRecords(20, 300, func(_, j int) byte { return byte('I' - j*30/300) })},
		// Low qualities ('!' and '#') exercise the quality-skip path fixed in stegano 0.6.
		{"low_quality", fastqRecords(60, 100, func(i, j int) byte { return "!#+5?I"[(i+j)%6] })},
	}

	e := newTestEngine(t)
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			fp, err := e.GenerateFingerprint("dataset-001")
			require.NoError(t, err)

			out, err := e.EmbedInQualities(tc.data, fp)
			require.NoError(t, err)

			got, err := e.ExtractFingerprint(out)
			require.NoError(t, err)
			assert.Equal(t, fp, got)
			assert.True(t, e.VerifyFingerprint(out))

			// Only quality scores may change, each by at most one Phred step.
			require.Equal(t, len(tc.data), len(out))
			in, mod := strings.Split(string(tc.data), "\n"), strings.Split(string(out), "\n")
			require.Equal(t, len(in), len(mod))
			changed := 0
			for i := range in {
				if i%4 != 3 {
					require.Equal(t, in[i], mod[i], "non-quality line %d changed", i)
					continue
				}
				for j := range in[i] {
					d := int(in[i][j]) - int(mod[i][j])
					require.True(t, d >= -1 && d <= 1, "line %d pos %d: quality %q -> %q", i, j, in[i][j], mod[i][j])
					if d != 0 {
						changed++
					}
				}
			}
			assert.Greater(t, changed, 0, "embedding changed no quality scores")
		})
	}
}

func TestSteganoEngine_CleanFASTQHasNoFingerprint(t *testing.T) {
	e := newTestEngine(t)
	clean := fastqRecords(50, 150, func(_, _ int) byte { return 'F' })

	_, err := e.ExtractFingerprint(clean)
	requireCode(t, err, ErrorCodeCorruption)
	assert.False(t, e.VerifyFingerprint(clean))
}

// Damaged quality scores must make extraction fail, never return a wrong fingerprint.
func TestSteganoEngine_CorruptedFingerprintIsRejected(t *testing.T) {
	e := newTestEngine(t)
	out, err := e.EmbedInQualities(fastqRecords(50, 150, func(_, _ int) byte { return 'F' }), "customer-123:dataset-001")
	require.NoError(t, err)

	lines := strings.Split(string(out), "\n")
	for i := 3; i < 12; i += 4 {
		q := []byte(lines[i])
		for j := range q {
			q[j] ^= 1
		}
		lines[i] = string(q)
	}

	got, err := e.ExtractFingerprint([]byte(strings.Join(lines, "\n")))
	require.Error(t, err, "extracted %q from corrupted data", got)
	requireCode(t, err, ErrorCodeCorruption)
}

func TestSteganoEngine_InsufficientCapacity(t *testing.T) {
	e := newTestEngine(t)
	_, err := e.EmbedInQualities(fastqRecords(1, 150, func(_, _ int) byte { return 'F' }), "customer-123:dataset-001")

	fe := requireCode(t, err, ErrorCodeCapacity)
	assert.ErrorIs(t, err, ErrCapacityExceeded)
	assert.True(t, fe.Recoverable)
}

// EmbedInQualities needs complete FASTQ records. A single read's bare quality
// string (what genomics.Engine currently passes per read) is rejected.
func TestSteganoEngine_RejectsBareQualityLine(t *testing.T) {
	_, err := newTestEngine(t).EmbedInQualities([]byte(strings.Repeat("I", 150)), "customer-123:dataset-001")
	requireCode(t, err, ErrorCodeValidation)
	assert.ErrorIs(t, err, ErrInvalidFormat)
}

func TestSteganoEngine_InvalidInputs(t *testing.T) {
	e := newTestEngine(t)
	data := fastqRecords(50, 150, func(_, _ int) byte { return 'F' })

	for name, tc := range map[string]struct {
		data []byte
		fp   string
	}{
		"empty_data":        {nil, "fp"},
		"empty_fingerprint": {data, ""},
		"control_chars":     {data, "fp\x00id"},
		"too_long":          {data, strings.Repeat("x", MaxFingerprintLength+1)},
	} {
		t.Run(name, func(t *testing.T) {
			_, err := e.EmbedInQualities(tc.data, tc.fp)
			requireCode(t, err, ErrorCodeValidation)
		})
	}

	_, err := e.ExtractFingerprint(nil)
	requireCode(t, err, ErrorCodeValidation)

	_, err = NewSteganoEngine("")
	requireCode(t, err, ErrorCodeValidation)
}

func samRecords(n int) []byte {
	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n")
	for i := 0; i < n; i++ {
		fmt.Fprintf(&b, "read_%d\t0\tchr1\t%d\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n", i, i*10+1)
	}
	return []byte(b.String())
}

func TestSteganoEngine_SAMAuxiliaryRoundTrip(t *testing.T) {
	e := newTestEngine(t)
	sam := samRecords(50)
	const fp = "customer-123:dataset-001"

	out, err := e.EmbedInSAMRecords(sam, fp)
	require.NoError(t, err)

	got, err := e.ExtractFromSAMRecords(out)
	require.NoError(t, err)
	assert.Equal(t, fp, got)

	// The auxiliary strategy only appends optional tags: headers and the 11
	// mandatory fields of every record are unchanged.
	in := strings.Split(strings.TrimSuffix(string(sam), "\n"), "\n")
	mod := strings.Split(strings.TrimSuffix(string(out), "\n"), "\n")
	require.Equal(t, len(in), len(mod))
	tagged := 0
	for i := range in {
		if strings.HasPrefix(in[i], "@") {
			assert.Equal(t, in[i], mod[i])
			continue
		}
		fields := strings.Split(mod[i], "\t")
		require.GreaterOrEqual(t, len(fields), 11)
		assert.Equal(t, in[i], strings.Join(fields[:11], "\t"), "record %d mandatory fields changed", i)
		for _, f := range fields[11:] {
			if strings.HasPrefix(f, "zs:Z:") {
				tagged++
			}
		}
	}
	assert.Greater(t, tagged, 0, "no zs:Z: tags written")
}

func TestSteganoEngine_CleanSAMHasNoFingerprint(t *testing.T) {
	e := newTestEngine(t)
	_, err := e.ExtractFromSAMRecords(samRecords(50))
	requireCode(t, err, ErrorCodeCorruption)

	_, err = e.EmbedInSAMRecords(nil, "fp")
	requireCode(t, err, ErrorCodeValidation)
	_, err = e.EmbedInSAMRecords(samRecords(50), "")
	requireCode(t, err, ErrorCodeValidation)
}
