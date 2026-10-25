package genomics

// Variants planted from a VCF (variants_path / --variants, #32).
//
// With SequencerConfig.VariantsPath set, the engine reads the VCF once when
// it is created and shares the variants, read-only, with every output
// stream: file output (Generate, GenerateFASTQFile, Regenerate) and
// streamed output (GenerateToWriter) apply them by the same rules:
//
//   - with reference_path, variants with CHROM/POS are applied in reference
//     coordinates and their REF is checked against the reference (see
//     engine_reference.go); strict_ref makes a mismatch an error before
//     anything is written;
//   - without a reference, positioned variants place reads in a region:
//     variants_region (--region) or, without it, the span of the positioned
//     variants plus one read length on each side (all on one chromosome);
//   - without positions, variants are applied position-blind.
//
// Every allele is a Bernoulli draw, r.Float64() < AF, on the stream the read
// (single-end) or fragment (paired-end: the batch's pair stream, so both
// mates and both mate files carry the same allele) is drawn from, so output
// stays a pure function of the seed, configuration and VCF, the same for any
// number of workers.
//
// File output records the planted alleles in a variants table per output
// unit, <name>.variants.tsv next to the R1 (or single-end) file, listed in
// the manifest with role "truth"; the manifest records the VCF by name, size
// and SHA-256 (never its path), and Regenerate refuses a run whose VCF it was
// not given (RegenOptions.Variants).

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

// LoadVariantsVCF reads the variants of a VCF file (plain, or gzip/bgzip
// compressed when the name ends in .gz): CHROM, POS, ID, REF, the first ALT
// of each record, and the allele frequency AF= from INFO (default 0.5, the
// heterozygous frequency, when absent or outside [0, 1]). Records with REF or
// ALT "." are skipped; an unparseable POS makes the variant position-blind.
func LoadVariantsVCF(path string) ([]VariantRecord, error) {
	f, err := os.Open(path) // #nosec G304 -- the operator's VCF (--variants / variants_path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("opening gzip reader for %s: %w", path, err)
		}
		defer gz.Close()
		r = gz
	}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	var records []VariantRecord
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}
		pos, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			pos = 0 // unparseable position: position-blind
		}
		ref, alt := fields[3], fields[4]
		if ref == "." || alt == "." || ref == "" || alt == "" {
			continue
		}
		if i := strings.IndexByte(alt, ','); i >= 0 {
			alt = alt[:i] // the first ALT of a multi-allelic site
		}
		freq := 0.5
		if len(fields) >= 8 {
			freq = ParseAF(fields[7], freq)
		}
		id := fields[2]
		if id == "." {
			id = ""
		}
		records = append(records, VariantRecord{ID: id, Ref: ref, Alt: alt, Freq: freq, Chromosome: fields[0], Position: pos})
	}
	return records, sc.Err()
}

// ParseAF returns the first AF=<value> of a VCF INFO column, or def when it
// is absent, unparseable or outside [0, 1].
func ParseAF(info string, def float64) float64 {
	for _, field := range strings.Split(info, ";") {
		val, ok := strings.CutPrefix(field, "AF=")
		if !ok {
			continue
		}
		if i := strings.IndexByte(val, ','); i >= 0 {
			val = val[:i]
		}
		if f, err := strconv.ParseFloat(val, 64); err == nil && f >= 0 && f <= 1 {
			return f
		}
	}
	return def
}

// ParseRegion parses chr:start-end (1-based, inclusive; commas allowed) into
// a 0-based half-open region.
func ParseRegion(s string) (*GenomicRegion, error) {
	bad := fmt.Errorf("invalid region %q: want chr:start-end (1-based, inclusive)", s)
	i := strings.LastIndexByte(s, ':')
	if i <= 0 {
		return nil, bad
	}
	startStr, endStr, ok := strings.Cut(strings.ReplaceAll(s[i+1:], ",", ""), "-")
	if !ok {
		return nil, bad
	}
	start, err1 := strconv.ParseInt(startStr, 10, 64)
	end, err2 := strconv.ParseInt(endStr, 10, 64)
	if err1 != nil || err2 != nil || start < 1 || end < start {
		return nil, bad
	}
	return &GenomicRegion{Chromosome: s[:i], Start: start - 1, End: end}, nil
}

// ResolveVariantRegion returns the region reads are placed in without a
// reference: region (chr:start-end) when given; otherwise, when variants
// carry positions, the span of the positioned variants extended by
// readLength on each side (clamped at 0); otherwise nil (position-blind
// injection). Positioned variants on more than one chromosome need region.
func ResolveVariantRegion(variants []VariantRecord, region string, readLength int) (*GenomicRegion, error) {
	if region != "" {
		r, err := ParseRegion(region)
		if err != nil {
			return nil, err
		}
		if r.End-r.Start <= int64(readLength) {
			return nil, fmt.Errorf("region %s is %d bases; it must be longer than the read length %d", region, r.End-r.Start, readLength)
		}
		return r, nil
	}
	var span *GenomicRegion
	chroms := map[string]bool{}
	for _, v := range variants {
		if v.Chromosome == "" || v.Position <= 0 {
			continue
		}
		chroms[v.Chromosome] = true
		start, end := v.Position-1, v.Position-1+int64(len(v.Ref))
		if span == nil {
			span = &GenomicRegion{Chromosome: v.Chromosome, Start: start, End: end}
			continue
		}
		span.Start = min(span.Start, start)
		span.End = max(span.End, end)
	}
	if len(chroms) > 1 {
		names := make([]string, 0, len(chroms))
		for c := range chroms {
			names = append(names, c)
		}
		sort.Strings(names)
		return nil, fmt.Errorf("variants are on %d chromosomes (%s): pass --region chr:start-end (variants_region) to choose where reads are placed", len(names), strings.Join(names, ", "))
	}
	if span != nil {
		span.Start = max(0, span.Start-int64(readLength))
		span.End += int64(readLength)
	}
	return span, nil
}

// loadVariants reads the configured VCF, records it for the manifest and
// prepares its variants: in reference mode checked against the reference
// (refVars), otherwise with the region reads are placed in. Errors are
// reported when generation starts (checkExtensions).
func (e *Engine) loadVariants() {
	c := e.config
	if c.VariantsRegion != "" && c.ReferencePath != "" {
		e.varErr = fmt.Errorf("variants_region %q places reads without a reference: with reference_path, set reference_region (--region) instead", c.VariantsRegion)
		return
	}
	if c.VariantsPath != "" {
		vs, err := LoadVariantsVCF(c.VariantsPath)
		if err != nil {
			e.varErr = fmt.Errorf("loading variants from %s: %w", c.VariantsPath, err)
			return
		}
		n, sum, err := manifest.HashFile(c.VariantsPath)
		if err != nil {
			e.varErr = fmt.Errorf("variants: %w", err)
			return
		}
		e.variants = vs
		e.variantsInput = &manifest.Input{Name: filepath.Base(c.VariantsPath), Bytes: n, SHA256: sum}
		fmt.Fprintf(e.logw, "   [info] Loaded %d variants from %s\n", len(vs), c.VariantsPath)
	}
	if e.ref != nil {
		if e.refVars, e.varErr = e.referenceVariants(e.variants); e.varErr == nil {
			e.refVars.report(e)
		}
		return
	}
	if e.refErr != nil {
		return
	}
	region, err := ResolveVariantRegion(e.variants, c.VariantsRegion, c.ReadLength)
	if err != nil {
		e.varErr = err
		return
	}
	e.varRegion = region
	if region != nil {
		e.reportUnplaced(region)
	}
}

// reportUnplaced logs the placement region and the variants region-mode
// injection cannot place: those without a position or on another
// chromosome.
func (e *Engine) reportUnplaced(region *GenomicRegion) {
	unpositioned, elsewhere := 0, 0
	for _, v := range e.variants {
		switch {
		case v.Chromosome == "" || v.Position <= 0:
			unpositioned++
		case !sameChromosome(v.Chromosome, region.Chromosome):
			elsewhere++
		}
	}
	fmt.Fprintf(e.logw, "   [info] Placing reads in %s:%d-%d\n", region.Chromosome, region.Start+1, region.End)
	if unpositioned > 0 {
		fmt.Fprintf(e.logw, "   [warn] %d variants without CHROM/POS are ignored in region mode\n", unpositioned)
	}
	if elsewhere > 0 {
		fmt.Fprintf(e.logw, "   [warn] %d variants on chromosomes other than %s are ignored\n", elsewhere, region.Chromosome)
	}
}

// Variants returns the variants loaded from variants_path (nil for none).
// The slice is shared: do not modify it.
func (e *Engine) Variants() []VariantRecord { return e.variants }

// VariantsInput returns the VCF the engine plants variants from, as the
// manifest records it, or nil.
func (e *Engine) VariantsInput() *manifest.Input { return e.variantsInput }

// configuredOpts returns the variant options of the engine's configuration
// (variants_path, variants_region): the variants and the placement region.
func (e *Engine) configuredOpts() GenerateOpts {
	return GenerateOpts{Variants: e.variants, RefRegion: e.varRegion}
}

// countOverlapSkips adds n skipped overlapping variant applications to the
// run's count (a flow cell's engine counts in its parent's).
func (e *Engine) countOverlapSkips(n int) {
	for e.parent != nil {
		e = e.parent
	}
	e.overlapSkips.Add(int64(n))
}

// plantedAllele is one variant allele written into a read or fragment: the
// variant, whether the allele is ALT, and its offset and length in the
// sequence (forward reference orientation).
type plantedAllele struct {
	v        *VariantRecord
	alt      bool
	off, len int
}

// plantLog collects the alleles written into one read (pair). A nil
// *plantLog records nothing.
type plantLog struct{ alleles []plantedAllele }

// splice records an allele of v written at off over refLen bases. Alleles
// already recorded to its right shift by the length change; one it
// overwrites is dropped.
func (l *plantLog) splice(v *VariantRecord, alt bool, off, refLen, alleleLen int) {
	if l == nil {
		return
	}
	kept := l.alleles[:0]
	for _, a := range l.alleles {
		switch {
		case a.off >= off+refLen:
			a.off += alleleLen - refLen
		case a.off+a.len > off:
			continue // overwritten (position-blind variants only)
		}
		kept = append(kept, a)
	}
	l.alleles = append(kept, plantedAllele{v: v, alt: alt, off: off, len: alleleLen})
}

// trim drops the alleles that start beyond the first n bases (an allele
// that starts inside them and runs past the end, such as an insertion at
// the end of a read, is kept: the read carries its first bases).
func (l *plantLog) trim(n int) {
	if l == nil {
		return
	}
	kept := l.alleles[:0]
	for _, a := range l.alleles {
		if a.off < n {
			kept = append(kept, a)
		}
	}
	l.alleles = kept
}

// remap moves each recorded allele to the offset f gives it (f(off, len)),
// dropping those f rejects: the alleles of a single-end read with indel
// errors move to their offsets in the read as written (#30).
func (l *plantLog) remap(f func(off, size int) (int, bool)) {
	if l == nil {
		return
	}
	kept := l.alleles[:0]
	for _, a := range l.alleles {
		if x, ok := f(a.off, a.len); ok {
			a.off = x
			kept = append(kept, a)
		}
	}
	l.alleles = kept
}

// reset clears the log for the next read (pair).
func (l *plantLog) reset() {
	if l != nil {
		l.alleles = l.alleles[:0]
	}
}

// variantsHeader is the first line of a variants table.
const variantsHeader = "#read_name\tread\tchrom\tpos\tid\tref\talt\tallele\toffset\n"

// writePlanted writes one variants-table line per allele planted in the
// read (pair) rec, sorted by offset: its name and 1-based number in the
// file (names repeat when reads are placed in a region without a
// reference), the variant (CHROM, 1-based POS, ID, REF, ALT; "." for none),
// which allele the read carries (REF or ALT), and the allele's 0-based
// offset in the fragment (paired-end) or read (single-end) on the
// reference's forward strand. An allele can run past the end of the read.
func writePlanted(w *bufio.Writer, rec *seqread.Record, l *plantLog) {
	sort.SliceStable(l.alleles, func(i, j int) bool { return l.alleles[i].off < l.alleles[j].off })
	dot := func(s string) string {
		if s == "" {
			return "."
		}
		return s
	}
	for _, a := range l.alleles {
		pos := "."
		if a.v.Position > 0 {
			pos = strconv.FormatInt(a.v.Position, 10)
		}
		allele := "REF"
		if a.alt {
			allele = "ALT"
		}
		// #nosec G104 -- bufio errors are sticky and surface at Flush
		w.WriteString(rec.Name + "\t" + strconv.Itoa(rec.Index) + "\t" + dot(a.v.Chromosome) + "\t" + pos + "\t" + dot(a.v.ID) + "\t" +
			a.v.Ref + "\t" + a.v.Alt + "\t" + allele + "\t" + strconv.Itoa(a.off) + "\n")
	}
}

// variantsPath returns the variants table written next to the FASTQ file
// fastq.
func variantsPath(fastq string) string {
	return strings.TrimSuffix(truthPath(fastq), ".truth.tsv") + ".variants.tsv"
}
