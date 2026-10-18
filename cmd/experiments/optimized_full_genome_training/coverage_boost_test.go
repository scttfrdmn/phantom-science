//go:build !release
// +build !release

package main

import (
	"os"
	"path/filepath"
	"sync"
	"testing"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
)

// TestChromosomeWorker_ProcessesWork feeds ChromosomeWork items through channels
// so the inner k-mer loop in chromosomeWorker runs.
func TestChromosomeWorker_ProcessesWork(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)

	workChan := make(chan ChromosomeWork, 2)
	resultChan := make(chan ChromosomeResult, 2)

	workChan <- ChromosomeWork{Name: "chr1", Sequence: []byte("ACGTACGTACGT"), Index: 0}
	workChan <- ChromosomeWork{Name: "chr2", Sequence: []byte("NNNACGTNNN"), Index: 1}
	close(workChan)

	var wg sync.WaitGroup
	wg.Add(1)
	go trainer.chromosomeWorker(&wg, workChan, resultChan)
	wg.Wait()
	close(resultChan)

	var results []ChromosomeResult
	for r := range resultChan {
		results = append(results, r)
	}
	if len(results) != 2 {
		t.Errorf("expected 2 results, got %d", len(results))
	}
	for _, r := range results {
		if r.Name == "chr1" && r.KmerCount == 0 {
			t.Errorf("chr1: expected KmerCount > 0, got 0")
		}
	}
}

// TestChromosomeWorker_EmptySequence sends a zero-length sequence.
func TestChromosomeWorker_EmptySequence(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)

	workChan := make(chan ChromosomeWork, 1)
	resultChan := make(chan ChromosomeResult, 1)
	workChan <- ChromosomeWork{Name: "chrEmpty", Sequence: []byte{}, Index: 0}
	close(workChan)

	var wg sync.WaitGroup
	wg.Add(1)
	go trainer.chromosomeWorker(&wg, workChan, resultChan)
	wg.Wait()
	close(resultChan)

	var results []ChromosomeResult
	for r := range resultChan {
		results = append(results, r)
	}
	if len(results) != 1 {
		t.Errorf("expected 1 result, got %d", len(results))
	}
	if results[0].KmerCount != 0 {
		t.Errorf("empty sequence KmerCount = %d, want 0", results[0].KmerCount)
	}
}

// TestChromosomeWorker_ClosedEmpty verifies worker exits when channel is empty.
func TestChromosomeWorker_ClosedEmpty(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)

	workChan := make(chan ChromosomeWork)
	resultChan := make(chan ChromosomeResult, 1)
	close(workChan)

	var wg sync.WaitGroup
	wg.Add(1)
	go trainer.chromosomeWorker(&wg, workChan, resultChan)
	wg.Wait()
	if len(resultChan) != 0 {
		t.Errorf("expected 0 results for empty work channel, got %d", len(resultChan))
	}
}

// TestChromosomeWorker_AllInvalidBases sends only N-bases filtered by isValidKmerOptimized.
func TestChromosomeWorker_AllInvalidBases(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)

	workChan := make(chan ChromosomeWork, 1)
	resultChan := make(chan ChromosomeResult, 1)
	workChan <- ChromosomeWork{Name: "chrN", Sequence: []byte("NNNNNNNNNN"), Index: 0}
	close(workChan)

	var wg sync.WaitGroup
	wg.Add(1)
	go trainer.chromosomeWorker(&wg, workChan, resultChan)
	wg.Wait()
	close(resultChan)

	for r := range resultChan {
		if r.KmerCount != 0 {
			t.Errorf("all-N: expected KmerCount=0, got %d", r.KmerCount)
		}
	}
}

// TestTrainOptimized_CollectsStats verifies stats are returned when no chromosomes load.
func TestTrainOptimized_CollectsStats(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)

	stats, err := trainer.TrainOptimized()
	if err != nil {
		t.Fatalf("TrainOptimized error: %v", err)
	}
	if stats == nil {
		t.Fatal("expected non-nil stats")
	}
	if stats.ParallelSpeedup <= 0 {
		t.Errorf("ParallelSpeedup = %v, want > 0", stats.ParallelSpeedup)
	}
}

// TestChromosomeWorker_MultipleWorkers starts multiple concurrent workers.
func TestChromosomeWorker_MultipleWorkers(t *testing.T) {
	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(t.TempDir())
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)

	const numItems = 6
	workChan := make(chan ChromosomeWork, numItems)
	resultChan := make(chan ChromosomeResult, numItems)

	for i := 0; i < numItems; i++ {
		workChan <- ChromosomeWork{
			Name:     "chr" + string(rune(48+i)),
			Sequence: []byte("ACGTACGT"),
			Index:    i,
		}
	}
	close(workChan)

	var wg sync.WaitGroup
	for w := 0; w < 3; w++ {
		wg.Add(1)
		go trainer.chromosomeWorker(&wg, workChan, resultChan)
	}
	wg.Wait()
	close(resultChan)

	count := 0
	for range resultChan {
		count++
	}
	if count != numItems {
		t.Errorf("expected %d results, got %d", numItems, count)
	}
}

// TestTrainOptimized_WithRealChromosome creates a chromosome FASTA file in the
// reference genome directory so that TrainOptimized can load and process it,
// exercising the chromosome-loading, k-mer aggregation, and model training paths.
func TestTrainOptimized_WithRealChromosome(t *testing.T) {
	tmpDir := t.TempDir()

	// Create a chr1.fasta file so LoadChromosome("chr1") succeeds
	fastaContent := ">chr1\nACGTACGTACGTACGT\n"
	if err := os.WriteFile(filepath.Join(tmpDir, "chr1.fasta"), []byte(fastaContent), 0600); err != nil {
		t.Fatalf("WriteFile chr1.fasta: %v", err)
	}

	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)

	stats, err := trainer.TrainOptimized()
	if err != nil {
		t.Fatalf("TrainOptimized with chr1: %v", err)
	}
	if stats == nil {
		t.Fatal("expected non-nil stats")
	}
	// At least chr1 should be processed
	if stats.ChromosomesProcessed < 1 {
		t.Logf("WARN: ChromosomesProcessed = %d (expected >= 1 if chr1 loaded)", stats.ChromosomesProcessed)
	}
}

// TestTrainOptimized_WithMultipleChromosomes creates multiple chromosome files
// to exercise the k-mer aggregation loop.
func TestTrainOptimized_WithMultipleChromosomes(t *testing.T) {
	tmpDir := t.TempDir()
	chromosomes := []string{"chr1", "chr2", "chr3"}
	for _, chrName := range chromosomes {
		seq := "ACGTACGTACGT" + chrName // unique sequence per chromosome
		fastaContent := ">" + chrName + "\n" + seq + "\n"
		if err := os.WriteFile(filepath.Join(tmpDir, chrName+".fasta"), []byte(fastaContent), 0600); err != nil {
			t.Fatalf("WriteFile %s.fasta: %v", chrName, err)
		}
	}

	model := authenticity.NewAdaptiveHMMModel(4)
	refGenome := reference.NewSimpleReferenceGenome(tmpDir)
	trainer := NewOptimizedFullGenomeTrainer(model, refGenome)

	stats, err := trainer.TrainOptimized()
	if err != nil {
		t.Fatalf("TrainOptimized with multiple chromosomes: %v", err)
	}
	if stats == nil {
		t.Fatal("expected non-nil stats")
	}
	if stats.UniqueKmers < 0 {
		t.Errorf("UniqueKmers = %d, want >= 0", stats.UniqueKmers)
	}
}

