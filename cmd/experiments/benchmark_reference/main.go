// Command benchmark_reference is a developer benchmark, not part of the
// phantom-science CLI. It times region reads and k-mer-profile lookups
// through reference.ReferenceManager, first against an empty cache and then
// re-reading the same regions, and reports per-phase time, errors and cache
// hits and misses (taken from the loader's own counters).
//
// Regions are tiled end to end along the chromosomes the reference's index
// lists, cycling through the chromosomes, so every region lies inside its
// chromosome whatever reference is used. Without -reference it writes a
// temporary FASTA of -gen-chroms chromosomes of -gen-length seeded random
// bases and deletes it afterwards.
//
// Benchmarks (-benchmark):
//
//	cache     read every region cold, then again warm
//	kmer      build every region's k-mer profile cold, then fetch it again
//	prefetch  start ReferenceManager.PrefetchRegions in the background, then read
//	all       cache, kmer, then prefetch, clearing the cache between them
//
// -prefetch turns on the manager's read-ahead, which widens every read by
// -prefetch-window bases each side. The loader caches by exact coordinates,
// so with read-ahead on, regions loaded by PrefetchRegions or WarmupReference
// are not the ones GetSequence asks for; the hit counts show this.
//
// Usage:
//
//	go run ./cmd/experiments/benchmark_reference -benchmark all -regions 500
//	go run ./cmd/experiments/benchmark_reference -reference hg38.fa -benchmark cache
package main

import (
	"bufio"
	"errors"
	"flag"
	"fmt"
	"io"
	"math/rand/v2"
	"os"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// maxK bounds -k: each profile is a dense 4^k packed table.
const maxK = 12

type config struct {
	reference      string
	cacheSize      int64
	k              int
	regions        int
	regionSize     int64
	benchmark      string
	verbose        bool
	prefetch       bool
	prefetchWindow int64
	genChroms      int
	genLength      int64
	seed           uint64
}

func main() {
	var cfg config
	flag.StringVar(&cfg.reference, "reference", "", "reference FASTA (default: generate a temporary one)")
	flag.Int64Var(&cfg.cacheSize, "cache", 2<<30, "cache size in bytes")
	flag.IntVar(&cfg.k, "k", 8, fmt.Sprintf("k-mer size for profiles (1-%d)", maxK))
	flag.IntVar(&cfg.regions, "regions", 1000, "number of regions per phase")
	flag.Int64Var(&cfg.regionSize, "region-size", 10000, "region length in bases")
	flag.StringVar(&cfg.benchmark, "benchmark", "all", "cache, kmer, prefetch or all")
	flag.BoolVar(&cfg.verbose, "verbose", false, "print the first error of each phase")
	flag.BoolVar(&cfg.prefetch, "prefetch", false, "enable the manager's read-ahead window")
	flag.Int64Var(&cfg.prefetchWindow, "prefetch-window", 10000, "read-ahead bases each side of a region when -prefetch is set")
	flag.IntVar(&cfg.genChroms, "gen-chroms", 5, "chromosomes in the generated reference")
	flag.Int64Var(&cfg.genLength, "gen-length", 10<<20, "bases per chromosome in the generated reference")
	flag.Uint64Var(&cfg.seed, "seed", 1, "seed for the generated reference")
	flag.Parse()

	if err := run(cfg, os.Stdout); err != nil {
		fmt.Fprintf(os.Stderr, "benchmark_reference: %v\n", err)
		os.Exit(1)
	}
}

func (c config) validate() error {
	switch {
	case c.k < 1 || c.k > maxK:
		return fmt.Errorf("-k must be between 1 and %d, got %d", maxK, c.k)
	case c.regions < 1:
		return fmt.Errorf("-regions must be at least 1, got %d", c.regions)
	case c.regionSize < 1:
		return fmt.Errorf("-region-size must be at least 1, got %d", c.regionSize)
	case c.cacheSize < 0:
		return fmt.Errorf("-cache must not be negative, got %d", c.cacheSize)
	case c.reference == "" && (c.genChroms < 1 || c.genLength < 1):
		return fmt.Errorf("-gen-chroms and -gen-length must be at least 1")
	}
	switch c.benchmark {
	case "cache", "kmer", "prefetch", "all":
	default:
		return fmt.Errorf("unknown -benchmark %q (want cache, kmer, prefetch or all)", c.benchmark)
	}
	return nil
}

func run(cfg config, w io.Writer) error {
	if err := cfg.validate(); err != nil {
		return err
	}

	refPath := cfg.reference
	if refPath == "" {
		p, err := createTestReference(cfg.genChroms, cfg.genLength, cfg.seed)
		if err != nil {
			return err
		}
		defer os.Remove(p)
		defer os.Remove(p + ".fai")
		refPath = p
		fmt.Fprintf(w, "generated reference: %s (%d x %d bases, seed %d)\n", p, cfg.genChroms, cfg.genLength, cfg.seed)
	}

	chroms, err := readLayout(refPath)
	if err != nil {
		return err
	}
	regions := planRegions(chroms, cfg.regions, cfg.regionSize)

	fmt.Fprintf(w, "reference: %s (%d chromosomes)\n", refPath, len(chroms))
	fmt.Fprintf(w, "regions: %d x %d bases, cache %.1f MiB, k=%d, read-ahead %t\n\n",
		len(regions), cfg.regionSize, float64(cfg.cacheSize)/(1024*1024), cfg.k, cfg.prefetch)

	manager, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: refPath,
		DefaultCacheSize:     cfg.cacheSize,
		EnablePrefetching:    cfg.prefetch,
		PrefetchWindow:       cfg.prefetchWindow,
		EnablePackedArrays:   true,
	})
	if err != nil {
		return fmt.Errorf("create reference manager: %w", err)
	}
	defer manager.Close()

	var phases []phase
	if cfg.benchmark == "cache" || cfg.benchmark == "all" {
		phases = append(phases, runCacheBenchmark(manager, regions)...)
	}
	if cfg.benchmark == "kmer" || cfg.benchmark == "all" {
		phases = append(phases, runKmerBenchmark(manager, regions, cfg.k)...)
	}
	if cfg.benchmark == "prefetch" || cfg.benchmark == "all" {
		phases = append(phases, runPrefetchBenchmark(manager, regions))
	}

	printPhases(w, phases, cfg.verbose)
	printCacheState(w, manager)

	for _, p := range phases {
		if p.Errors == p.Requests {
			return fmt.Errorf("every request in phase %q failed: %v", p.Name, p.FirstErr)
		}
	}
	return nil
}

// chromInfo is one chromosome of the reference.
type chromInfo struct {
	Name   string
	Length int64
}

// readLayout returns the chromosomes of the FASTA at path, in file order,
// from its .fai index or a scan when there is none.
func readLayout(path string) ([]chromInfo, error) {
	parser, err := reference.NewFastaParser(path)
	if err != nil {
		return nil, fmt.Errorf("open reference %s: %w", path, err)
	}
	index, err := parser.LoadIndex()
	if err != nil {
		return nil, fmt.Errorf("index reference %s: %w", path, err)
	}
	var chroms []chromInfo
	for _, name := range index.OrderedNames {
		if info := index.Chromosomes[name]; info != nil && info.Length > 0 {
			chroms = append(chroms, chromInfo{Name: name, Length: info.Length})
		}
	}
	if len(chroms) == 0 {
		return nil, fmt.Errorf("reference %s has no non-empty sequences", path)
	}
	return chroms, nil
}

// planRegions returns n regions of size bases. Region i is on chromosome
// i mod len(chroms); successive regions on a chromosome are adjacent tiles,
// wrapping to the start when they reach the end. A chromosome shorter than
// size yields regions covering the whole chromosome.
func planRegions(chroms []chromInfo, n int, size int64) []reference.GenomicRegion {
	regions := make([]reference.GenomicRegion, n)
	for i := range regions {
		c := chroms[i%len(chroms)]
		tile := int64(i / len(chroms))
		var start int64
		if tiles := c.Length / size; tiles > 0 {
			start = (tile % tiles) * size
		}
		end := min(start+size, c.Length)
		regions[i] = reference.GenomicRegion{Chromosome: c.Name, Start: start, End: end}
	}
	return regions
}

// createTestReference writes a FASTA of chroms chromosomes (chr1, chr2, ...)
// of length seeded random bases to a temporary file and returns its path.
func createTestReference(chroms int, length int64, seed uint64) (string, error) {
	f, err := os.CreateTemp("", "benchmark_reference_*.fasta")
	if err != nil {
		return "", fmt.Errorf("create temporary reference: %w", err)
	}
	werr := writeTestReference(f, chroms, length, seed)
	cerr := f.Close()
	if err := errors.Join(werr, cerr); err != nil {
		os.Remove(f.Name())
		return "", fmt.Errorf("write temporary reference: %w", err)
	}
	return f.Name(), nil
}

// writeTestReference writes chroms FASTA records of length uniformly random
// bases, 80 per line.
func writeTestReference(w io.Writer, chroms int, length int64, seed uint64) error {
	const bases, lineLen = "ACGT", 80
	rng := rand.New(rand.NewPCG(seed, 0x9e3779b97f4a7c15)) // #nosec G404 -- benchmark input, not security-sensitive
	bw := bufio.NewWriterSize(w, 1<<20)
	line := make([]byte, lineLen+1)
	for c := 1; c <= chroms; c++ {
		fmt.Fprintf(bw, ">chr%d\n", c)
		for done := int64(0); done < length; {
			n := int(min(lineLen, length-done))
			for i := 0; i < n; i++ {
				line[i] = bases[rng.IntN(4)]
			}
			line[n] = '\n'
			if _, err := bw.Write(line[:n+1]); err != nil {
				return err
			}
			done += int64(n)
		}
	}
	return bw.Flush()
}

// phase is one timed pass over the regions.
type phase struct {
	Name     string
	Requests int
	Errors   int
	FirstErr error
	Duration time.Duration
	Hits     int64 // loader cache hits during the pass
	Misses   int64 // loader cache misses during the pass
	Kmers    int64 // k-mers in the returned profiles (kmer phases)
}

func loaderCounters(m *reference.ReferenceManager) (hits, misses int64) {
	s := m.GetCombinedStatistics()
	return s.TotalCacheHits, s.TotalCacheMisses
}

// timePass runs fn for each region and records time, errors and the change
// in the loader's hit and miss counters.
func timePass(m *reference.ReferenceManager, name string, regions []reference.GenomicRegion,
	fn func(reference.GenomicRegion) (int64, error)) phase {
	p := phase{Name: name, Requests: len(regions)}
	h0, m0 := loaderCounters(m)
	start := time.Now()
	for _, r := range regions {
		n, err := fn(r)
		if err != nil {
			p.Errors++
			if p.FirstErr == nil {
				p.FirstErr = fmt.Errorf("%s:%d-%d: %w", r.Chromosome, r.Start, r.End, err)
			}
			continue
		}
		p.Kmers += n
	}
	p.Duration = time.Since(start)
	h1, m1 := loaderCounters(m)
	p.Hits, p.Misses = h1-h0, m1-m0
	return p
}

func readSequence(m *reference.ReferenceManager) func(reference.GenomicRegion) (int64, error) {
	return func(r reference.GenomicRegion) (int64, error) {
		_, err := m.GetSequence(r.Chromosome, r.Start, r.End)
		return 0, err
	}
}

// runCacheBenchmark reads every region into an empty cache, then reads them
// all again.
func runCacheBenchmark(m *reference.ReferenceManager, regions []reference.GenomicRegion) []phase {
	m.ClearAllCaches()
	cold := timePass(m, "cache: cold read", regions, readSequence(m))
	warm := timePass(m, "cache: warm read", regions, readSequence(m))
	return []phase{cold, warm}
}

// runKmerBenchmark builds every region's k-mer profile from an empty cache,
// then requests each profile again.
func runKmerBenchmark(m *reference.ReferenceManager, regions []reference.GenomicRegion, k int) []phase {
	m.ClearAllCaches()
	profile := func(r reference.GenomicRegion) (int64, error) {
		p, err := m.GetKmerProfile(r.Chromosome, r.Start, r.End, k)
		if err != nil {
			return 0, err
		}
		return int64(p.TotalCount), nil
	}
	cold := timePass(m, "kmer: build profile", regions, profile)
	warm := timePass(m, "kmer: cached profile", regions, profile)
	return []phase{cold, warm}
}

// runPrefetchBenchmark starts PrefetchRegions on an empty cache and reads the
// regions while it runs; hits and misses include the background loader's.
// The background loader is not waited for, so this phase runs last.
func runPrefetchBenchmark(m *reference.ReferenceManager, regions []reference.GenomicRegion) phase {
	m.ClearAllCaches()
	started := false
	read := readSequence(m)
	return timePass(m, "prefetch: read during PrefetchRegions", regions, func(r reference.GenomicRegion) (int64, error) {
		if !started {
			started = true
			if err := m.PrefetchRegions(regions, 0); err != nil {
				return 0, fmt.Errorf("PrefetchRegions: %w", err)
			}
		}
		return read(r)
	})
}

func printPhases(w io.Writer, phases []phase, verbose bool) {
	fmt.Fprintf(w, "%-38s %8s %7s %12s %12s %8s %8s\n", "phase", "regions", "errors", "time", "regions/s", "hits", "misses")
	for _, p := range phases {
		rate := "n/a"
		if s := p.Duration.Seconds(); s > 0 {
			rate = fmt.Sprintf("%.0f", float64(p.Requests-p.Errors)/s)
		}
		fmt.Fprintf(w, "%-38s %8d %7d %12s %12s %8d %8d\n", p.Name, p.Requests, p.Errors,
			p.Duration.Round(time.Microsecond), rate, p.Hits, p.Misses)
		if p.Kmers > 0 {
			fmt.Fprintf(w, "%-38s k-mers counted: %d\n", "", p.Kmers)
		}
		if verbose && p.FirstErr != nil {
			fmt.Fprintf(w, "%-38s first error: %v\n", "", p.FirstErr)
		}
	}
	fmt.Fprintln(w)
}

func printCacheState(w io.Writer, m *reference.ReferenceManager) {
	s, ok := m.GetCombinedStatistics().LoaderStats["default"]
	if !ok {
		return
	}
	fmt.Fprintf(w, "cache after last phase: %d regions, %.1f MiB of %.1f MiB (%.1f%%)\n",
		s.CachedRegions, float64(s.CacheSize)/(1024*1024), float64(s.MaxCacheSize)/(1024*1024),
		s.CacheUtilization*100)
}
