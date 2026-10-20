package fingerprinting

import (
	"hash/maphash"
	"math"
)

// Read matching for keyed-seed verification.
//
// A Matcher indexes the reads of suspect files, then is shown every read of
// a regenerated run (Observe). A suspect read matches when its sequence, or
// its reverse complement, equals a regenerated read (exact), or when at
// least half of up to MatchSlots k-mers sampled along it (MatchK bases, either
// strand) occur in regenerated reads (contained). Containment finds reads
// that were trimmed, re-oriented or had a few bases edited; qualities and
// headers are never used, so rewritten qualities and renamed reads still
// match.
//
// Memory is proportional to the suspect reads, not to the regenerated run.

// MatchK is the k-mer length used for containment matching.
const MatchK = 31

// MatchSlots is the most k-mers sampled from one suspect read.
const MatchSlots = 32

const (
	kmerMask   = 1<<(2*MatchK) - 1
	filterBits = 24
)

// slotRef names one sampled k-mer of one suspect read.
type slotRef struct {
	read uint32
	slot uint8
}

// suspectRead is one indexed suspect read.
type suspectRead struct {
	file   int32
	length int32
	slots  uint8  // sampled k-mers
	found  uint32 // bit per slot seen in the regenerated run
	exact  bool
}

// FileMatch is the outcome for one suspect file.
type FileMatch struct {
	Name string
	// Checked counts the suspect reads examined.
	Checked int
	// Exact counts reads equal to a regenerated read (either strand).
	Exact int
	// Contained counts reads that are not exact but share at least half
	// of their sampled k-mers with the regenerated run.
	Contained int
	// ExpectedChance bounds the matches expected by chance if the
	// regenerated and suspect reads were independent uniform-random
	// sequences. It does not hold for low-complexity sequence or reads
	// sampled from a shared reference genome; measure those with a decoy
	// run instead.
	ExpectedChance float64
}

// Matched is the number of matching reads.
func (f FileMatch) Matched() int { return f.Exact + f.Contained }

// Matcher matches suspect reads against a regenerated run.
type Matcher struct {
	seed   maphash.Seed
	files  []string
	reads  []suspectRead
	exact  map[uint64][]uint32
	kmers  map[uint64][]slotRef
	filter []uint64
	rc     []byte

	// Regenerated reads, bases and k-mer positions observed since Reset.
	regenReads, regenBases, regenKmers int64
}

// NewMatcher returns an empty Matcher.
func NewMatcher() *Matcher {
	return &Matcher{
		seed:   maphash.MakeSeed(),
		exact:  map[uint64][]uint32{},
		kmers:  map[uint64][]slotRef{},
		filter: make([]uint64, 1<<filterBits/64),
	}
}

// AddFile registers a suspect file and returns its index for AddRead.
func (m *Matcher) AddFile(name string) int {
	m.files = append(m.files, name)
	return len(m.files) - 1
}

// AddRead indexes one suspect read of file. Lower-case bases are folded to
// upper case; bases other than ACGT break k-mers.
func (m *Matcher) AddRead(file int, seq []byte) {
	s := make([]byte, len(seq))
	for i, b := range seq {
		if 'a' <= b && b <= 'z' {
			b -= 'a' - 'A'
		}
		s[i] = b
	}
	id := uint32(len(m.reads))                                 // #nosec G115 -- suspect reads are capped far below 2^32
	r := suspectRead{file: int32(file), length: int32(len(s))} // #nosec G115
	h := m.canonicalHash(s)
	m.exact[h] = append(m.exact[h], id)

	if n := len(s) - MatchK + 1; n > 0 {
		slots := min(n, MatchSlots)
		for i := 0; i < slots; i++ {
			pos := 0
			if slots > 1 {
				pos = i * (n - 1) / (slots - 1)
			}
			km, ok := canonicalKmer(s[pos : pos+MatchK])
			if !ok {
				continue
			}
			m.kmers[km] = append(m.kmers[km], slotRef{read: id, slot: r.slots})
			m.setFilter(km)
			r.slots++
		}
	}
	m.reads = append(m.reads, r)
}

// Reads returns the number of indexed suspect reads.
func (m *Matcher) Reads() int { return len(m.reads) }

// Reset clears the matches, keeping the index, for another regenerated run.
func (m *Matcher) Reset() {
	for i := range m.reads {
		m.reads[i].found, m.reads[i].exact = 0, false
	}
	m.regenReads, m.regenBases, m.regenKmers = 0, 0, 0
}

// Observe matches one regenerated read against the index.
func (m *Matcher) Observe(seq []byte) {
	m.regenReads++
	m.regenBases += int64(len(seq))
	if ids, ok := m.exact[m.canonicalHash(seq)]; ok {
		for _, id := range ids {
			if int(m.reads[id].length) == len(seq) {
				m.reads[id].exact = true
			}
		}
	}
	var fwd, rev uint64
	valid := 0
	for _, b := range seq {
		c, ok := baseCode(b)
		if !ok {
			valid = 0
			continue
		}
		fwd = (fwd<<2 | c) & kmerMask
		rev = rev>>2 | (3-c)<<(2*(MatchK-1))
		if valid++; valid < MatchK {
			continue
		}
		m.regenKmers++
		km := min(fwd, rev)
		if !m.inFilter(km) {
			continue
		}
		for _, ref := range m.kmers[km] {
			m.reads[ref.read].found |= 1 << ref.slot
		}
	}
}

// Regenerated returns the reads and bases observed since Reset.
func (m *Matcher) Regenerated() (reads, bases int64) { return m.regenReads, m.regenBases }

// Results returns one FileMatch per registered file, in AddFile order.
func (m *Matcher) Results() []FileMatch {
	out := make([]FileMatch, len(m.files))
	for i, name := range m.files {
		out[i].Name = name
	}
	// p is the chance that one sampled k-mer occurs among the regenerated
	// k-mers (either strand) of uniform-random sequence.
	p := math.Min(1, 2*float64(m.regenKmers)/math.Pow(4, MatchK))
	for _, r := range m.reads {
		f := &out[r.file]
		f.Checked++
		switch {
		case r.exact:
			f.Exact++
		case r.slots > 0 && 2*popcount(r.found) >= int(r.slots):
			f.Contained++
		}
		f.ExpectedChance += chanceMatch(int(r.length), int(r.slots), p, m.regenReads)
	}
	return out
}

// chanceMatch bounds the probability that a suspect read of length n with
// slots sampled k-mers matches uniform-random regenerated sequence: an exact
// match among reads regenerated reads, or ceil(slots/2) sampled k-mers each
// present with probability p (union bound over the choice of k-mers).
func chanceMatch(n, slots int, p float64, reads int64) float64 {
	exact := math.Min(1, 2*float64(reads)*math.Pow(4, -float64(n)))
	if slots == 0 || p == 0 {
		return exact
	}
	need := (slots + 1) / 2
	lchoose, _ := math.Lgamma(float64(slots + 1))
	a, _ := math.Lgamma(float64(need + 1))
	b, _ := math.Lgamma(float64(slots - need + 1))
	contained := math.Exp(lchoose - a - b + float64(need)*math.Log(p))
	return math.Min(1, exact+contained)
}

func popcount(x uint32) int {
	n := 0
	for ; x != 0; x &= x - 1 {
		n++
	}
	return n
}

// canonicalHash hashes the lesser of seq and its reverse complement.
func (m *Matcher) canonicalHash(seq []byte) uint64 {
	if cap(m.rc) < len(seq) {
		m.rc = make([]byte, len(seq))
	}
	rc := m.rc[:len(seq)]
	for i, b := range seq {
		rc[len(seq)-1-i] = complement(b)
	}
	if string(rc) < string(seq) {
		return maphash.Bytes(m.seed, rc)
	}
	return maphash.Bytes(m.seed, seq)
}

func complement(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'T':
		return 'A'
	}
	return b
}

func baseCode(b byte) (uint64, bool) {
	switch b {
	case 'A':
		return 0, true
	case 'C':
		return 1, true
	case 'G':
		return 2, true
	case 'T':
		return 3, true
	}
	return 0, false
}

// canonicalKmer encodes kmer (MatchK bases) as the lesser of its strand
// encodings; ok is false if it has a base other than ACGT.
func canonicalKmer(kmer []byte) (uint64, bool) {
	var fwd, rev uint64
	for _, b := range kmer {
		c, ok := baseCode(b)
		if !ok {
			return 0, false
		}
		fwd = fwd<<2 | c
		rev = rev>>2 | (3-c)<<(2*(MatchK-1))
	}
	return min(fwd, rev), true
}

func filterIndex(km uint64) uint64 { return (km * 0x9e3779b97f4a7c15) >> (64 - filterBits) }

func (m *Matcher) setFilter(km uint64) {
	i := filterIndex(km)
	m.filter[i/64] |= 1 << (i % 64)
}

func (m *Matcher) inFilter(km uint64) bool {
	i := filterIndex(km)
	return m.filter[i/64]&(1<<(i%64)) != 0
}
