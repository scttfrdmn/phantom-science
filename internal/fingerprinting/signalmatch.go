package fingerprinting

import (
	"encoding/binary"
	"hash/maphash"
)

// Raw-signal matching for keyed-seed verification of BLOW5 output (#31).
//
// A SignalMatcher indexes the raw signal (int16 ADC samples) of suspect
// reads, then is shown the raw signal of every read of a regenerated signal
// run (Observe). A suspect read matches when its whole signal equals a
// regenerated read's (exact), or when at least half of up to SignalSlots
// windows of SignalWindow consecutive samples, sampled along it, occur in
// regenerated signal (contained). Read IDs, calibration fields and auxiliary
// fields are never used, so signal whose reads were renamed, re-packed
// (BLOW5, SLOW5, record or signal compression) or trimmed still matches.
//
// The samples of a read are a deterministic function of the run seed, the
// read ID and its sequence, so a window of SignalWindow samples carries the
// seed's noise; how much depends on the profile's noise and dwell model, so
// no fixed chance bound is given. Verification measures the background
// with a decoy run of the same inputs under another seed instead.
//
// Memory is proportional to the suspect reads, not to the regenerated run.

// SignalWindow is the number of consecutive samples in a matched window.
const SignalWindow = 16

// SignalSlots is the most windows sampled from one suspect read.
const SignalSlots = 32

// signalPrime is the multiplier of the polynomial window hash (odd).
const signalPrime = 0x100000001b3

// suspectSignal is one indexed suspect read.
type suspectSignal struct {
	file    int32
	samples int32
	slots   uint8
	found   uint32
	exact   bool
}

// SignalMatcher matches suspect raw signal against a regenerated run.
type SignalMatcher struct {
	seed    maphash.Seed
	files   []string
	reads   []suspectSignal
	exact   map[uint64][]uint32
	windows map[uint64][]slotRef
	filter  []uint64
	pow     uint64 // signalPrime^SignalWindow
	buf     []byte

	regenReads, regenSamples int64
}

// NewSignalMatcher returns an empty SignalMatcher.
func NewSignalMatcher() *SignalMatcher {
	pow := uint64(1)
	for range SignalWindow {
		pow *= signalPrime
	}
	return &SignalMatcher{
		seed:    maphash.MakeSeed(),
		exact:   map[uint64][]uint32{},
		windows: map[uint64][]slotRef{},
		filter:  make([]uint64, 1<<filterBits/64),
		pow:     pow,
	}
}

// AddFile registers a suspect file and returns its index for AddRead.
func (m *SignalMatcher) AddFile(name string) int {
	m.files = append(m.files, name)
	return len(m.files) - 1
}

// AddRead indexes the raw signal of one suspect read of file.
func (m *SignalMatcher) AddRead(file int, raw []int16) {
	id := uint32(len(m.reads))                                      // #nosec G115 -- suspect reads are capped far below 2^32
	r := suspectSignal{file: int32(file), samples: int32(len(raw))} // #nosec G115
	h := m.fullHash(raw)
	m.exact[h] = append(m.exact[h], id)
	if n := len(raw) - SignalWindow + 1; n > 0 {
		slots := min(n, SignalSlots)
		for i := 0; i < slots; i++ {
			pos := 0
			if slots > 1 {
				pos = i * (n - 1) / (slots - 1)
			}
			w := windowKey(windowHash(raw[pos : pos+SignalWindow]))
			m.windows[w] = append(m.windows[w], slotRef{read: id, slot: r.slots})
			m.setFilter(w)
			r.slots++
		}
	}
	m.reads = append(m.reads, r)
}

// Reads returns the number of indexed suspect reads.
func (m *SignalMatcher) Reads() int { return len(m.reads) }

// Reset clears the matches, keeping the index, for another regenerated run.
func (m *SignalMatcher) Reset() {
	for i := range m.reads {
		m.reads[i].found, m.reads[i].exact = 0, false
	}
	m.regenReads, m.regenSamples = 0, 0
}

// Observe matches the raw signal of one regenerated read against the index.
func (m *SignalMatcher) Observe(raw []int16) {
	m.regenReads++
	m.regenSamples += int64(len(raw))
	if ids, ok := m.exact[m.fullHash(raw)]; ok {
		for _, id := range ids {
			if int(m.reads[id].samples) == len(raw) {
				m.reads[id].exact = true
			}
		}
	}
	if len(raw) < SignalWindow {
		return
	}
	h := windowHash(raw[:SignalWindow])
	for i := SignalWindow; ; i++ {
		if w := windowKey(h); m.inFilter(w) {
			for _, ref := range m.windows[w] {
				m.reads[ref.read].found |= 1 << ref.slot
			}
		}
		if i == len(raw) {
			break
		}
		h = h*signalPrime + sample(raw[i]) - sample(raw[i-SignalWindow])*m.pow
	}
}

// Regenerated returns the reads and samples observed since Reset.
func (m *SignalMatcher) Regenerated() (reads, samples int64) { return m.regenReads, m.regenSamples }

// Results returns one FileMatch per registered file, in AddFile order.
// ExpectedChance is zero: the background is measured by a decoy run.
func (m *SignalMatcher) Results() []FileMatch {
	out := make([]FileMatch, len(m.files))
	for i, name := range m.files {
		out[i].Name = name
	}
	for _, r := range m.reads {
		f := &out[r.file]
		f.Checked++
		switch {
		case r.exact:
			f.Exact++
		case r.slots > 0 && 2*popcount(r.found) >= int(r.slots):
			f.Contained++
			if popcount(r.found) == int(r.slots) {
				f.Whole++
			}
		}
	}
	return out
}

func sample(v int16) uint64 { return uint64(uint16(v)) } // #nosec G115 -- bit pattern

// windowHash is the polynomial hash Σ s_j · P^(W-1-j) of a window, mod 2^64.
func windowHash(w []int16) uint64 {
	var h uint64
	for _, v := range w {
		h = h*signalPrime + sample(v)
	}
	return h
}

// windowKey scrambles a window hash (splitmix64 finalizer) so that its high
// bits, which index the filter, depend on every sample.
func windowKey(h uint64) uint64 {
	h ^= h >> 30
	h *= 0xbf58476d1ce4e5b9
	h ^= h >> 27
	h *= 0x94d049bb133111eb
	return h ^ h>>31
}

// fullHash hashes a whole signal.
func (m *SignalMatcher) fullHash(raw []int16) uint64 {
	m.buf = m.buf[:0]
	for _, v := range raw {
		m.buf = binary.LittleEndian.AppendUint16(m.buf, uint16(v)) // #nosec G115 -- bit pattern
	}
	return maphash.Bytes(m.seed, m.buf)
}

func (m *SignalMatcher) setFilter(w uint64) {
	i := filterIndex(w)
	m.filter[i/64] |= 1 << (i % 64)
}

func (m *SignalMatcher) inFilter(w uint64) bool {
	i := filterIndex(w)
	return m.filter[i/64]&(1<<(i%64)) != 0
}
