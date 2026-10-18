// Package nanopore models Oxford Nanopore (MinKNOW) run metadata for
// synthetic reads: run information, UUID4 read IDs, FASTQ header fields, a
// per-channel timing model, and the sequencing_summary.txt writer.
//
// It imports only the standard library and the seqread leaf package, so the
// genomics engine can use it without an import cycle.
package nanopore

import (
	"fmt"
	"math/rand/v2"
	"strings"
	"time"
	"unicode"
)

// Defaults for a PromethION R10.4.1 run.
const (
	// DefaultChannels is the number of channels on a PromethION flow cell
	// (FLO-PRO114M: "up to 2675 nanopore channels",
	// https://store.nanoporetech.com/flow-cells.html). MinION/GridION flow
	// cells have 512.
	DefaultChannels = 2675
	// DefaultSpeedBPS is the translocation speed of Kit 14 / R10.4.1
	// chemistry (400 bases per second, the MinKNOW default since 2023).
	DefaultSpeedBPS = 400.0
	// DefaultProtocolGroupID is the experiment name used when none is set.
	DefaultProtocolGroupID = "phantom_science"
	// DefaultSampleID is MinKNOW's sample ID when none is entered.
	DefaultSampleID = "no_sample"
)

// RunInfo is the run-level metadata written to every read header and to the
// sequencing summary.
type RunInfo struct {
	RunID           string    // 40 lowercase hex characters, like MinKNOW
	FlowCellID      string    // e.g. PAW12345
	ProtocolGroupID string    // MinKNOW experiment name
	SampleID        string    // MinKNOW sample ID ("no_sample" when unset, as MinKNOW)
	Start           time.Time // acquisition start; read start times are offsets from it
	Channels        int       // flow-cell channels; reads are numbered 1..Channels
	SpeedBPS        float64   // translocation speed, bases per second
}

// runInfoStream is the PCG stream used to derive default run metadata from
// the seed; it is independent of every read stream.
const runInfoStream = 0x6f6e745f72756e // "ont_run"

// WithDefaults returns ri with every unset field filled in. Generated values
// (run ID, flow-cell ID, start time) are derived only from seed, so a run is
// reproducible from its seed.
func (ri RunInfo) WithDefaults(seed uint64) RunInfo {
	// #nosec G404 -- simulation metadata, reproducible by design
	r := rand.New(rand.NewPCG(seed, runInfoStream))
	// Always draw in the same order so each default does not depend on
	// which other fields were set.
	runID := randomHex(r, 20)
	flowCell := fmt.Sprintf("PA%c%05d", 'A'+rune(r.IntN(26)), r.IntN(100000))
	start := time.Date(2024, 1, 1, 0, 0, 0, 0, time.UTC).Add(time.Duration(r.Int64N(365*24*3600)) * time.Second)

	if ri.RunID == "" {
		ri.RunID = runID
	}
	if ri.FlowCellID == "" {
		ri.FlowCellID = flowCell
	}
	if ri.SampleID == "" {
		ri.SampleID = DefaultSampleID
	}
	if ri.ProtocolGroupID == "" {
		ri.ProtocolGroupID = DefaultProtocolGroupID
	}
	if ri.Start.IsZero() {
		ri.Start = start
	}
	if ri.Channels == 0 {
		ri.Channels = DefaultChannels
	}
	if ri.SpeedBPS == 0 {
		ri.SpeedBPS = DefaultSpeedBPS
	}
	return ri
}

// Validate reports run metadata that would produce malformed headers or an
// unusable timing model. Header values must be single tokens.
func (ri RunInfo) Validate() error {
	for _, f := range []struct{ name, v string }{
		{"run_id", ri.RunID}, {"flow_cell_id", ri.FlowCellID},
		{"protocol_group_id", ri.ProtocolGroupID}, {"sample_id", ri.SampleID},
	} {
		if f.v == "" {
			return fmt.Errorf("ONT %s is empty", f.name)
		}
		if strings.IndexFunc(f.v, func(c rune) bool { return unicode.IsSpace(c) || c == '=' || !unicode.IsPrint(c) }) >= 0 {
			return fmt.Errorf("ONT %s %q must not contain whitespace or '='", f.name, f.v)
		}
	}
	if ri.Channels <= 0 {
		return fmt.Errorf("ONT channels must be positive, got %d", ri.Channels)
	}
	if !(ri.SpeedBPS > 0) {
		return fmt.Errorf("ONT speed_bps must be positive, got %g", ri.SpeedBPS)
	}
	return nil
}

// ShortRunID returns the first 8 characters of the run ID, as MinKNOW uses
// in file names.
func (ri RunInfo) ShortRunID() string {
	if len(ri.RunID) > 8 {
		return ri.RunID[:8]
	}
	return ri.RunID
}

const hexDigits = "0123456789abcdef"

// randomHex returns 2n lowercase hex characters drawn from r.
func randomHex(r *rand.Rand, n int) string {
	b := make([]byte, 0, 2*n)
	for i := 0; i < n; i++ {
		v := byte(r.Uint32())
		b = append(b, hexDigits[v>>4], hexDigits[v&15])
	}
	return string(b)
}

// NewUUID4 returns a random (version 4, RFC 4122 variant) UUID drawn from r,
// in the canonical 8-4-4-4-12 lowercase form. MinKNOW read IDs are UUID4s.
func NewUUID4(r *rand.Rand) string {
	var u [16]byte
	a, b := r.Uint64(), r.Uint64()
	for i := 0; i < 8; i++ {
		u[i] = byte(a >> (8 * i))
		u[8+i] = byte(b >> (8 * i))
	}
	u[6] = u[6]&0x0f | 0x40 // version 4
	u[8] = u[8]&0x3f | 0x80 // variant 10xx (RFC 4122)
	var s [36]byte
	j := 0
	for i, c := range u {
		if i == 4 || i == 6 || i == 8 || i == 10 {
			s[j] = '-'
			j++
		}
		s[j], s[j+1] = hexDigits[c>>4], hexDigits[c&15]
		j += 2
	}
	return string(s[:])
}
