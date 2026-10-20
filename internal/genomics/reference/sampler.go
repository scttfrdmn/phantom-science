package reference

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"math/rand/v2"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Contig is one sequence of a reference FASTA.
type Contig struct {
	Name   string
	Length int64
}

// maxFragmentAttempts bounds the position redraws Fragment makes inside the
// chosen contig before it gives up on a window free of N.
const maxFragmentAttempts = 10000

// Sampler draws reference fragments from an uncompressed FASTA file.
//
// NewSampler reads the sequence it needs into memory once: every contig, or
// only the requested region. Bases are stored upper-cased at one byte per
// base, so a whole GRCh38 primary assembly takes about 3.1 GB and a region
// takes its own length in bytes. Positions that are not A, C, G or T (N and
// the other IUPAC ambiguity codes) are also recorded as a sorted list of
// runs, a few bytes per run.
//
// After NewSampler returns, a Sampler is never modified, so its methods are
// safe for concurrent use without locks. Fragment draws from the *rand.Rand
// it is given; a rand.Rand is not safe for concurrent use, so give each
// goroutine its own.
//
// Only uncompressed FASTA is read. A samtools .fai index next to the file
// (path + ".fai") is used when present; otherwise the index is built by
// scanning the file. Gzip and bgzip files are rejected.
type Sampler struct {
	path    string
	contigs []Contig       // every sequence in the file, in file order
	byName  map[string]int // contig name -> index into contigs

	loaded   []loadedSeq
	loadedBy map[string]int // contig name -> index into loaded

	cum      []int64 // cumulative A/C/G/T counts over loaded, for contig choice
	minClean int64   // shortest maxClean over loaded
}

// loadedSeq is a contig, or the requested part of one, held in memory.
type loadedSeq struct {
	name   string
	offset int64  // 0-based contig coordinate of seq[0]
	seq    []byte // upper-cased bases
	acgt   int64  // number of A/C/G/T bases in seq
	// maxClean is the longest run of A/C/G/T, the longest fragment the
	// sequence can supply.
	maxClean int64
	// Runs of bases that are not A/C/G/T, as [badStart[i], badEnd[i]) in
	// seq coordinates, sorted and disjoint.
	badStart []int64
	badEnd   []int64
}

var (
	upperTable [256]byte
	isLetter   [256]bool
	isACGT     [256]bool
	complement [256]byte
)

func init() {
	for i := 0; i < 256; i++ {
		c := byte(i)
		if c >= 'a' && c <= 'z' {
			c -= 'a' - 'A'
		}
		upperTable[i] = c
		isLetter[i] = c >= 'A' && c <= 'Z'
	}
	for _, c := range []byte("ACGT") {
		isACGT[c] = true
	}
	complement['A'], complement['C'], complement['G'], complement['T'] = 'T', 'G', 'C', 'A'
}

// NewSampler indexes the FASTA at path and loads its sequence into memory.
//
// region is "" for the whole file, a contig name ("chr1"), or
// "contig:start-end" with 1-based inclusive coordinates; commas are allowed
// ("chr1:1,000,001-2,000,000") and "contig:start" runs to the contig end.
// Contig names are matched with or without a "chr" prefix (chr1 and 1, chrM
// and MT). With a region, only that part of the contig is loaded and
// Fragment and Bases are confined to it.
func NewSampler(path, region string) (*Sampler, error) {
	if err := rejectCompressed(path); err != nil {
		return nil, err
	}
	index, err := loadOrBuildIndex(path)
	if err != nil {
		return nil, err
	}
	if len(index.OrderedNames) == 0 {
		return nil, fmt.Errorf("reference %s contains no sequences", path)
	}

	s := &Sampler{
		path:     path,
		byName:   make(map[string]int, len(index.OrderedNames)),
		loadedBy: make(map[string]int),
	}
	for _, name := range index.OrderedNames {
		info := index.Chromosomes[name]
		s.byName[name] = len(s.contigs)
		s.contigs = append(s.contigs, Contig{Name: name, Length: info.Length})
	}

	f, err := os.Open(path) // #nosec G304 -- path is the user's reference
	if err != nil {
		return nil, fmt.Errorf("open reference: %w", err)
	}
	defer f.Close()

	if region == "" {
		for _, name := range index.OrderedNames {
			info := index.Chromosomes[name]
			if info.Length == 0 {
				continue
			}
			if err := s.load(f, info, 0, info.Length); err != nil {
				return nil, err
			}
		}
	} else {
		name, start, end, err := s.parseRegion(region)
		if err != nil {
			return nil, err
		}
		if err := s.load(f, index.Chromosomes[name], start, end); err != nil {
			return nil, err
		}
	}
	if len(s.loaded) == 0 {
		return nil, fmt.Errorf("reference %s contains no bases", path)
	}

	var total int64
	s.minClean = -1
	for _, l := range s.loaded {
		total += l.acgt
		s.cum = append(s.cum, total)
		if s.minClean < 0 || l.maxClean < s.minClean {
			s.minClean = l.maxClean
		}
	}
	if total == 0 {
		return nil, fmt.Errorf("reference %s has no A/C/G/T bases in the requested sequence", path)
	}
	return s, nil
}

// Contigs returns every sequence in the FASTA, in file order, with its full
// length, whether or not a region restricts sampling (this is the list a
// SAM header's @SQ lines need).
func (s *Sampler) Contigs() []Contig {
	return append([]Contig(nil), s.contigs...)
}

// Fragment draws one n-base fragment.
//
// The draw order is fixed, so the same *rand.Rand state and the same FASTA
// (and region) always give the same fragment:
//
//  1. contig: one r.Int64N(W), where W is the total A/C/G/T count of the
//     eligible loaded sequences, those with a run of at least n bases free
//     of N; each eligible contig is chosen with probability proportional
//     to its non-N (A/C/G/T) length.
//  2. position: one r.Int64N(L-n+1) per attempt, L being the loaded length
//     of that contig; a window containing a base other than A/C/G/T is
//     redrawn in the same contig, up to 10000 attempts, then an error is
//     returned. The accepted position is uniform over the contig's clean
//     windows.
//  3. strand: one r.IntN(2); 1 is the minus strand.
//
// chrom is the FASTA's own contig name and start the 0-based forward-strand
// coordinate of the fragment's leftmost base, whatever the strand. On the
// minus strand seq is the reverse complement of the forward bases
// [start, start+n). seq is upper-case A/C/G/T and owned by the caller.
func (s *Sampler) Fragment(r *rand.Rand, n int) (chrom string, start int64, minus bool, seq []byte, err error) {
	if n <= 0 {
		return "", 0, false, nil, fmt.Errorf("fragment length must be positive, got %d", n)
	}
	n64 := int64(n)

	cum := s.cum
	idx := []int(nil) // maps a cum position to s.loaded when some are excluded
	if n64 > s.minClean {
		var total int64
		cum = nil
		for i, l := range s.loaded {
			if l.maxClean >= n64 {
				total += l.acgt
				cum = append(cum, total)
				idx = append(idx, i)
			}
		}
		if len(cum) == 0 {
			return "", 0, false, nil, fmt.Errorf("fragment length %d is longer than every N-free stretch of %s (longest %d bases)", n, s.describeLoaded(), s.maxCleanLen())
		}
	}

	u := r.Int64N(cum[len(cum)-1])
	k := sort.Search(len(cum), func(i int) bool { return cum[i] > u })
	if idx != nil {
		k = idx[k]
	}
	l := &s.loaded[k]

	windows := int64(len(l.seq)) - n64 + 1
	pos := int64(-1)
	for attempt := 0; attempt < maxFragmentAttempts; attempt++ {
		p := r.Int64N(windows)
		if l.clean(p, p+n64) {
			pos = p
			break
		}
	}
	if pos < 0 {
		return "", 0, false, nil, fmt.Errorf("no %d-base window free of N found in %s after %d attempts", n, l.name, maxFragmentAttempts)
	}
	minus = r.IntN(2) == 1

	src := l.seq[pos : pos+n64]
	seq = make([]byte, n)
	if minus {
		for i, c := range src {
			seq[n-1-i] = complement[c]
		}
	} else {
		copy(seq, src)
	}
	return l.name, l.offset + pos, minus, seq, nil
}

// Bases returns the upper-cased forward-strand bases [start, end) (0-based,
// half-open) of chrom, which may be given with or without a "chr" prefix.
// The range must lie inside the loaded sequence: the region, when NewSampler
// was given one. N and other IUPAC codes are returned as they are.
func (s *Sampler) Bases(chrom string, start, end int64) ([]byte, error) {
	name, ok := s.resolve(chrom)
	if !ok {
		return nil, fmt.Errorf("contig %q not found in %s", chrom, s.path)
	}
	k, ok := s.loadedBy[name]
	if !ok {
		return nil, fmt.Errorf("contig %s is outside the loaded region (%s)", name, s.describeLoaded())
	}
	l := &s.loaded[k]
	lo, hi := l.offset, l.offset+int64(len(l.seq))
	if start < lo || end > hi || start >= end {
		return nil, fmt.Errorf("invalid range %s:%d-%d (loaded %s:%d-%d, 0-based half-open)", name, start, end, name, lo, hi)
	}
	return append([]byte(nil), l.seq[start-l.offset:end-l.offset]...), nil
}

// clean reports whether [a, b) of the loaded sequence is all A/C/G/T.
func (l *loadedSeq) clean(a, b int64) bool {
	i := sort.Search(len(l.badEnd), func(i int) bool { return l.badEnd[i] > a })
	return i == len(l.badEnd) || l.badStart[i] >= b
}

func (s *Sampler) maxCleanLen() int64 {
	var m int64
	for _, l := range s.loaded {
		m = max(m, l.maxClean)
	}
	return m
}

func (s *Sampler) describeLoaded() string {
	if len(s.loaded) == 1 {
		l := s.loaded[0]
		if l.offset != 0 || int64(len(l.seq)) != s.contigs[s.byName[l.name]].Length {
			return fmt.Sprintf("region %s:%d-%d", l.name, l.offset+1, l.offset+int64(len(l.seq)))
		}
	}
	return s.path
}

// resolve maps a contig name, with or without a "chr" prefix, to the name
// used in the FASTA.
func (s *Sampler) resolve(name string) (string, bool) {
	candidates := []string{name}
	switch {
	case name == "chrM":
		candidates = append(candidates, "MT", "M")
	case name == "MT" || name == "M":
		candidates = append(candidates, "chrM")
	case strings.HasPrefix(name, "chr"):
		candidates = append(candidates, name[3:])
	default:
		candidates = append(candidates, "chr"+name)
	}
	for _, c := range candidates {
		if _, ok := s.byName[c]; ok {
			return c, true
		}
	}
	return "", false
}

// parseRegion turns a 1-based inclusive region string into a FASTA contig
// name and 0-based half-open coordinates.
func (s *Sampler) parseRegion(region string) (name string, start, end int64, err error) {
	if n, ok := s.resolve(region); ok { // a whole contig (names may contain ':')
		return n, 0, s.contigs[s.byName[n]].Length, nil
	}
	colon := strings.LastIndexByte(region, ':')
	if colon < 0 {
		return "", 0, 0, fmt.Errorf("region %q: contig not found in %s", region, s.path)
	}
	n, ok := s.resolve(region[:colon])
	if !ok {
		return "", 0, 0, fmt.Errorf("region %q: contig %q not found in %s", region, region[:colon], s.path)
	}
	length := s.contigs[s.byName[n]].Length
	coords := strings.ReplaceAll(region[colon+1:], ",", "")
	from, to, hasDash := strings.Cut(coords, "-")
	first, err := strconv.ParseInt(from, 10, 64)
	if err != nil || first < 1 {
		return "", 0, 0, fmt.Errorf("region %q: start must be a positive integer (1-based)", region)
	}
	last := length
	if hasDash && to != "" {
		last, err = strconv.ParseInt(to, 10, 64)
		if err != nil {
			return "", 0, 0, fmt.Errorf("region %q: end must be an integer", region)
		}
	}
	if last < first {
		return "", 0, 0, fmt.Errorf("region %q: end is before start", region)
	}
	if last > length {
		return "", 0, 0, fmt.Errorf("region %q: end %d is past the end of %s (length %d)", region, last, n, length)
	}
	return n, first - 1, last, nil
}

// load reads bases [start, end) of info's contig into memory.
func (s *Sampler) load(f *os.File, info *ChromosomeInfo, start, end int64) error {
	seq, err := readBases(f, info, start, end)
	if err != nil {
		return err
	}
	l := loadedSeq{name: info.Name, offset: start, seq: seq}
	inRun := false
	for i, c := range seq {
		switch {
		case isACGT[c]:
			l.acgt++
			if inRun {
				l.badEnd = append(l.badEnd, int64(i))
				inRun = false
			}
		case !inRun:
			l.badStart = append(l.badStart, int64(i))
			inRun = true
		}
	}
	if inRun {
		l.badEnd = append(l.badEnd, int64(len(seq)))
	}
	prev := int64(0) // end of the previous bad run
	for i := range l.badStart {
		l.maxClean = max(l.maxClean, l.badStart[i]-prev)
		prev = l.badEnd[i]
	}
	l.maxClean = max(l.maxClean, int64(len(seq))-prev)
	s.loadedBy[info.Name] = len(s.loaded)
	s.loaded = append(s.loaded, l)
	return nil
}

// readBases reads bases [start, end) of a contig from the FASTA via its
// index, dropping line terminators and upper-casing. Anything else that is
// not a letter (a header, say, because the .fai does not match the file) is
// an error.
func readBases(f *os.File, info *ChromosomeInfo, start, end int64) ([]byte, error) {
	ll, bpl := int64(info.LineLength), int64(info.BytesPerLine)
	offset := func(p int64) int64 { return info.StartOffset + (p/ll)*bpl + p%ll }
	from, to := offset(start), offset(end-1)+1

	seq := make([]byte, end-start)
	got := 0
	buf := make([]byte, 4<<20)
	for pos := from; pos < to; {
		chunk := buf[:min(int64(len(buf)), to-pos)]
		m, err := f.ReadAt(chunk, pos)
		if m < len(chunk) {
			if err == nil || errors.Is(err, io.EOF) {
				err = io.ErrUnexpectedEOF
			}
			return nil, fmt.Errorf("read %s from %s: %w (is the .fai index stale?)", info.Name, f.Name(), err)
		}
		for _, c := range chunk {
			if c == '\n' || c == '\r' {
				continue
			}
			if !isLetter[c] || got == len(seq) {
				return nil, fmt.Errorf("%s: unexpected byte %q in %s at file offset %d; the .fai index does not match the FASTA (re-run samtools faidx)", f.Name(), c, info.Name, pos)
			}
			seq[got] = upperTable[c]
			got++
		}
		pos += int64(len(chunk))
	}
	if got != len(seq) {
		return nil, fmt.Errorf("%s: read %d of %d bases of %s; the .fai index does not match the FASTA (re-run samtools faidx)", f.Name(), got, len(seq), info.Name)
	}
	return seq, nil
}

// rejectCompressed fails for gzip and bgzip files, which the sampler does
// not read.
func rejectCompressed(path string) error {
	f, err := os.Open(path) // #nosec G304 -- path is the user's reference
	if err != nil {
		return fmt.Errorf("open reference: %w", err)
	}
	defer f.Close()
	var magic [2]byte
	if n, _ := io.ReadFull(f, magic[:]); n == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
		return fmt.Errorf("reference %s is gzip/bgzip-compressed; only uncompressed FASTA is supported: decompress it (gunzip, or bgzip -d) and optionally index it with samtools faidx", path)
	}
	return nil
}

// loadOrBuildIndex parses path+".fai" when it exists and otherwise indexes
// the FASTA by scanning it.
func loadOrBuildIndex(path string) (*FastaIndex, error) {
	fai := path + ".fai"
	if _, err := os.Stat(fai); err != nil {
		if !errors.Is(err, os.ErrNotExist) {
			return nil, fmt.Errorf("stat %s: %w", fai, err)
		}
		return scanFastaIndex(path)
	}
	parser := &FastaParser{filename: path, indexFile: fai}
	index, err := parser.LoadIndex()
	if err != nil {
		return nil, err
	}
	for _, name := range index.OrderedNames {
		info := index.Chromosomes[name]
		if info.Length > 0 && (info.LineLength <= 0 || info.BytesPerLine < info.LineLength) {
			return nil, fmt.Errorf("%s: malformed entry for %s", fai, name)
		}
	}
	if len(index.OrderedNames) != len(index.Chromosomes) {
		return nil, fmt.Errorf("%s: duplicate sequence names", fai)
	}
	return index, nil
}

// scanFastaIndex builds a samtools-compatible index by reading the FASTA
// once. Byte offsets are exact (CRLF line endings and lines of any length
// are handled). As with samtools faidx, every line of a sequence but the
// last must have the same length, since random access depends on it.
func scanFastaIndex(path string) (*FastaIndex, error) {
	f, err := os.Open(path) // #nosec G304 -- path is the user's reference
	if err != nil {
		return nil, fmt.Errorf("open reference: %w", err)
	}
	defer f.Close()

	index := &FastaIndex{Chromosomes: make(map[string]*ChromosomeInfo)}
	br := bufio.NewReaderSize(f, 1<<20)
	var (
		offset int64
		cur    *ChromosomeInfo
		ended  bool // cur has had a short (final) line
	)
	for {
		var (
			raw, content int64
			first        byte
			header       []byte
			last2        [2]byte
			eof          bool
		)
		for {
			chunk, err := br.ReadSlice('\n')
			if raw == 0 && len(chunk) > 0 {
				first = chunk[0]
			}
			if first == '>' && len(header) < 64<<10 {
				header = append(header, chunk...)
			}
			for _, c := range chunk {
				last2[0], last2[1] = last2[1], c
			}
			raw += int64(len(chunk))
			if err == bufio.ErrBufferFull {
				continue
			}
			if errors.Is(err, io.EOF) {
				eof = true
				break
			}
			if err != nil {
				return nil, fmt.Errorf("read %s: %w", path, err)
			}
			break
		}
		if raw == 0 && eof {
			break
		}
		content = raw
		if last2[1] == '\n' {
			content--
			if raw >= 2 && last2[0] == '\r' {
				content--
			}
		} else if last2[1] == '\r' {
			content--
		}
		lineStart := offset
		offset += raw

		switch {
		case first == '>':
			fields := strings.Fields(string(header[1:]))
			if len(fields) == 0 {
				return nil, fmt.Errorf("%s: header with no sequence name at byte %d", path, lineStart)
			}
			name := fields[0]
			if _, dup := index.Chromosomes[name]; dup {
				return nil, fmt.Errorf("%s: duplicate sequence name %s", path, name)
			}
			cur = &ChromosomeInfo{Name: name, StartOffset: offset}
			index.Chromosomes[name] = cur
			index.OrderedNames = append(index.OrderedNames, name)
			ended = false
		case content == 0:
			if cur != nil && cur.Length > 0 {
				ended = true
			} else if cur != nil {
				cur.StartOffset = offset // blank line before the first sequence line
			}
		case cur == nil:
			return nil, fmt.Errorf("%s: sequence before the first header", path)
		default:
			if cur.Length == 0 {
				cur.LineLength, cur.BytesPerLine = int(content), int(raw)
			} else if ended || content > int64(cur.LineLength) {
				return nil, fmt.Errorf("%s: sequence %s has lines of different lengths; rewrap it (e.g. seqkit seq -w 60) so random access works", path, cur.Name)
			}
			if content < int64(cur.LineLength) || (raw != int64(cur.BytesPerLine) && !eof) {
				ended = true
			}
			cur.Length += content
		}
		if eof {
			break
		}
	}
	return index, nil
}
