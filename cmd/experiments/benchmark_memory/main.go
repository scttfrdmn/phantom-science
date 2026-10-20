// Command benchmark_memory is a developer benchmark, not part of the
// phantom-science CLI. It counts k-mers in seeded random sequences with the
// hash-map counter (kmer.FrequencyExtractor) and the packed-array counter
// (kmer.PackedFrequencyExtractor) and reports, for each, the wall time,
// throughput, heap retained by the result and bytes allocated while
// counting, then checks that the two give identical counts.
//
// Retained heap is the change in runtime.MemStats.HeapAlloc between a GC
// before the run and a GC after it, with the result still live; allocated
// bytes is the change in TotalAlloc. Both are process-wide, so other
// goroutines' allocations are included.
//
// The packed counter allocates a dense 4^k table per input sequence, so k is
// limited to 1-12 (a 64 MiB table per sequence at k=12).
//
// Usage:
//
//	go run ./cmd/experiments/benchmark_memory -k 8 -sequences 5000 -length 1000
package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"math/rand/v2"
	"os"
	"runtime"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/kmer"
)

// maxK bounds k: the packed counter allocates 4^k 32-bit counters per
// sequence before merging.
const maxK = 12

// memoryScanKs are the k values the -memory table covers.
var memoryScanKs = []int{4, 6, 8, 10}

type config struct {
	k          int
	sequences  int
	length     int
	seed       uint64
	verbose    bool
	packedOnly bool
	mapOnly    bool
	memory     bool
	accuracy   bool
}

func main() {
	var cfg config
	flag.IntVar(&cfg.k, "k", 8, fmt.Sprintf("k-mer size (1-%d)", maxK))
	flag.IntVar(&cfg.sequences, "sequences", 1000, "number of generated sequences")
	flag.IntVar(&cfg.length, "length", 1000, "length of each generated sequence")
	flag.Uint64Var(&cfg.seed, "seed", 1, "seed for the generated sequences")
	flag.BoolVar(&cfg.verbose, "verbose", false, "print table statistics")
	flag.BoolVar(&cfg.packedOnly, "packed-only", false, "run the packed-array counter only")
	flag.BoolVar(&cfg.mapOnly, "map-only", false, "run the hash-map counter only")
	flag.BoolVar(&cfg.memory, "memory", false, fmt.Sprintf("also print retained heap for k in %v", memoryScanKs))
	flag.BoolVar(&cfg.accuracy, "accuracy", false, "print the five most frequent k-mers from both counters")
	flag.Parse()

	if err := run(cfg, os.Stdout); err != nil {
		fmt.Fprintf(os.Stderr, "benchmark_memory: %v\n", err)
		os.Exit(1)
	}
}

func (c config) validate() error {
	switch {
	case c.k < 1 || c.k > maxK:
		return fmt.Errorf("-k must be between 1 and %d, got %d", maxK, c.k)
	case c.sequences < 1:
		return fmt.Errorf("-sequences must be at least 1, got %d", c.sequences)
	case c.length < c.k:
		return fmt.Errorf("-length (%d) must be at least -k (%d)", c.length, c.k)
	case c.packedOnly && c.mapOnly:
		return errors.New("-packed-only and -map-only are mutually exclusive")
	case c.accuracy && (c.packedOnly || c.mapOnly):
		return errors.New("-accuracy needs both counters; drop -packed-only/-map-only")
	}
	return nil
}

func run(cfg config, w io.Writer) error {
	if err := cfg.validate(); err != nil {
		return err
	}

	fmt.Fprintf(w, "k-mer counter memory benchmark\n")
	fmt.Fprintf(w, "sequences=%d length=%d k=%d seed=%d (%.2f MiB of sequence)\n\n",
		cfg.sequences, cfg.length, cfg.k, cfg.seed,
		float64(cfg.sequences*cfg.length)/(1024*1024))

	seqs := generateTestSequences(cfg.sequences, cfg.length, cfg.seed)
	totalBases := cfg.sequences * cfg.length

	var (
		mapTable    *kmer.FrequencyTable
		packedTable *kmer.PackedFrequencyTable
		mapMeas     measurement
		packedMeas  measurement
	)

	if !cfg.packedOnly {
		var err error
		mapTable, mapMeas, err = countMap(seqs, cfg.k)
		if err != nil {
			return fmt.Errorf("hash-map counter: %w", err)
		}
		fmt.Fprintf(w, "Hash map\n")
		printMeasurement(w, mapMeas, totalBases, mapTable.TotalCount, len(mapTable.Frequencies))
		if cfg.verbose {
			s := mapTable.GetKmerStats()
			fmt.Fprintf(w, "  complexity=%.4f shannon=%.4f gc=%.4f mean=%.2f max=%d\n",
				s.Complexity, s.ShannonEntropy, s.GCContent, s.MeanCount, s.MaxCount)
		}
		fmt.Fprintln(w)
	}

	if !cfg.mapOnly {
		var err error
		packedTable, packedMeas, err = countPacked(seqs, cfg.k)
		if err != nil {
			return fmt.Errorf("packed-array counter: %w", err)
		}
		fmt.Fprintf(w, "Packed array\n")
		printMeasurement(w, packedMeas, totalBases, packedTable.TotalCount, packedTable.NonZeroCount)
		if cfg.verbose {
			s := packedTable.GetStatistics()
			fmt.Fprintf(w, "  table=%.2f MiB bits/item=%d words=%d max=%d sparsity=%.4f\n",
				float64(s.MemoryUsage)/(1024*1024), s.BitsPerItem, s.PackedArraySize,
				s.MaxFrequency, s.Sparsity)
		}
		fmt.Fprintln(w)
	}

	if mapTable != nil && packedTable != nil {
		fmt.Fprintf(w, "Comparison\n")
		fmt.Fprintf(w, "  time map/packed:     %s\n", ratio(float64(mapMeas.Duration), float64(packedMeas.Duration)))
		fmt.Fprintf(w, "  retained map/packed: %s\n", ratio(float64(mapMeas.RetainedBytes), float64(packedMeas.RetainedBytes)))
		verr := verifyResults(mapTable, packedTable, cfg.k)
		fmt.Fprintf(w, "  results match: %t\n", verr == nil)
		if verr != nil {
			fmt.Fprintf(w, "  mismatch: %v\n", verr)
		}
		fmt.Fprintln(w)

		if cfg.accuracy {
			fmt.Fprintf(w, "Top 5 k-mers (map count / packed count)\n")
			for i, kf := range mapTable.GetTopKmers(5) {
				pc, _ := packedTable.GetFrequency(kf.Kmer)
				fmt.Fprintf(w, "  %d. %s %d / %d\n", i+1, kf.Kmer, kf.Count, pc)
			}
			fmt.Fprintln(w)
		}
		if verr != nil {
			return fmt.Errorf("counters disagree: %w", verr)
		}
	}

	if cfg.memory {
		if err := runMemoryScan(w, seqs); err != nil {
			return err
		}
	}
	return nil
}

// generateTestSequences returns count sequences of length uniformly random
// bases from a PCG generator seeded with seed, so runs are reproducible.
func generateTestSequences(count, length int, seed uint64) [][]byte {
	const bases = "ACGT"
	rng := rand.New(rand.NewPCG(seed, 0x9e3779b97f4a7c15)) // #nosec G404 -- benchmark input, not security-sensitive
	seqs := make([][]byte, count)
	for i := range seqs {
		s := make([]byte, length)
		for j := range s {
			s[j] = bases[rng.IntN(4)]
		}
		seqs[i] = s
	}
	return seqs
}

// measurement is one counter run.
type measurement struct {
	Duration       time.Duration
	RetainedBytes  int64  // HeapAlloc after GC, minus before; may be negative
	AllocatedBytes uint64 // TotalAlloc delta
}

// measure runs fn between two forced GCs and reports its wall time, the
// heap still held once it returns (with its result kept live), and the
// bytes it allocated.
func measure[T any](fn func() (T, error)) (T, measurement, error) {
	var before, after runtime.MemStats
	runtime.GC()
	runtime.ReadMemStats(&before)

	start := time.Now()
	v, err := fn()
	elapsed := time.Since(start)

	runtime.GC()
	runtime.ReadMemStats(&after)
	runtime.KeepAlive(v)

	return v, measurement{
		Duration:       elapsed,
		RetainedBytes:  int64(after.HeapAlloc) - int64(before.HeapAlloc), // #nosec G115 -- heap sizes fit in int64
		AllocatedBytes: after.TotalAlloc - before.TotalAlloc,
	}, err
}

func countMap(seqs [][]byte, k int) (*kmer.FrequencyTable, measurement, error) {
	return measure(func() (*kmer.FrequencyTable, error) {
		return kmer.NewFrequencyExtractor(k).ExtractFromMultipleSequences(seqs)
	})
}

func countPacked(seqs [][]byte, k int) (*kmer.PackedFrequencyTable, measurement, error) {
	return measure(func() (*kmer.PackedFrequencyTable, error) {
		return kmer.NewPackedFrequencyExtractor(k).ExtractFromMultipleSequences(seqs)
	})
}

func printMeasurement(w io.Writer, m measurement, totalBases, totalKmers, uniqueKmers int) {
	secs := m.Duration.Seconds()
	fmt.Fprintf(w, "  duration:     %v\n", m.Duration)
	if secs > 0 {
		fmt.Fprintf(w, "  throughput:   %.2f MiB/s, %.0f k-mers/s\n",
			float64(totalBases)/(1024*1024)/secs, float64(totalKmers)/secs)
	}
	fmt.Fprintf(w, "  k-mers:       %d total, %d distinct\n", totalKmers, uniqueKmers)
	fmt.Fprintf(w, "  retained:     %.2f MiB\n", float64(m.RetainedBytes)/(1024*1024))
	fmt.Fprintf(w, "  allocated:    %.2f MiB\n", float64(m.AllocatedBytes)/(1024*1024))
}

// ratio formats a/b, or "n/a" when b is not positive (a retained-heap delta
// can be zero or negative when GC reclaims more than the run kept).
func ratio(a, b float64) string {
	if b <= 0 || a < 0 {
		return "n/a"
	}
	return fmt.Sprintf("%.2fx", a/b)
}

func runMemoryScan(w io.Writer, seqs [][]byte) error {
	fmt.Fprintf(w, "Retained heap by k\n")
	fmt.Fprintf(w, "%-4s %14s %14s %10s\n", "k", "map MiB", "packed MiB", "map/packed")
	for _, k := range memoryScanKs {
		if len(seqs) > 0 && len(seqs[0]) < k {
			continue
		}
		_, mm, err := countMap(seqs, k)
		if err != nil {
			return fmt.Errorf("hash-map counter, k=%d: %w", k, err)
		}
		_, pm, err := countPacked(seqs, k)
		if err != nil {
			return fmt.Errorf("packed-array counter, k=%d: %w", k, err)
		}
		fmt.Fprintf(w, "%-4d %14.2f %14.2f %10s\n", k,
			float64(mm.RetainedBytes)/(1024*1024), float64(pm.RetainedBytes)/(1024*1024),
			ratio(float64(mm.RetainedBytes), float64(pm.RetainedBytes)))
	}
	fmt.Fprintln(w)
	return nil
}

// verifyResults reports the first difference between the two tables, or nil
// when they hold the same total and the same count for every k-mer.
func verifyResults(mapResult *kmer.FrequencyTable, packedResult *kmer.PackedFrequencyTable, k int) error {
	if mapResult.TotalCount != packedResult.TotalCount {
		return fmt.Errorf("total count: map %d, packed %d", mapResult.TotalCount, packedResult.TotalCount)
	}
	if len(mapResult.Frequencies) != packedResult.NonZeroCount {
		return fmt.Errorf("distinct k-mers: map %d, packed %d", len(mapResult.Frequencies), packedResult.NonZeroCount)
	}
	for s, mc := range mapResult.Frequencies {
		pc, err := packedResult.GetFrequency(s)
		if err != nil {
			return fmt.Errorf("k-mer %s: packed lookup: %w", s, err)
		}
		if pc != mc {
			return fmt.Errorf("k-mer %s: map %d, packed %d", s, mc, pc)
		}
	}
	for i := 0; i < packedResult.MaxKmers; i++ {
		pc := packedResult.PackedFrequencies.Get(i)
		if pc == 0 {
			continue
		}
		s := kmer.IndexToKmerString(i, k)
		if mc, ok := mapResult.Frequencies[s]; !ok || mc != int(pc) {
			return fmt.Errorf("k-mer %s: map %d, packed %d", s, mc, pc)
		}
	}
	return nil
}
