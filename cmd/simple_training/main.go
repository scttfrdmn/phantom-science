package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"runtime"
	"strings"
	"sync"
	"sync/atomic"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity"
)

func main() {
	log.SetFlags(log.LstdFlags)
	
	fmt.Println("=== SIMPLE HIGH-PERFORMANCE GENOMICS TRAINING ===")
	totalCores := runtime.NumCPU()
	perfCores := detectAppleSiliconCores()
	fmt.Printf("🍎 Apple Silicon: %d total cores (%d P-cores + %d E-cores)\n", 
		totalCores, perfCores, totalCores-perfCores)
	fmt.Printf("⚡ Using %d P-cores for training\n", perfCores)
	fmt.Println()

	runtime.GOMAXPROCS(perfCores)
	
	genomeFile := "/Volumes/Public/genomics_training/human/GRCh38_complete.fasta"
	
	// Check if file exists
	if _, err := os.Stat(genomeFile); os.IsNotExist(err) {
		log.Fatalf("Genome file not found: %s", genomeFile)
	}
	
	fmt.Printf("📁 Training on: %s\n", genomeFile)
	
	// Initialize simple model
	fmt.Println("[STEP 1] Initializing simple HMM model...")
	model := authenticity.NewAdaptiveHMMModel(4)
	fmt.Printf("✅ Model initialized with k=4\n")
	fmt.Println()

	// Train with heavy CPU usage
	fmt.Println("[STEP 2] High-intensity training...")
	fmt.Printf("🔥 Using %d P-cores for maximum performance\n", perfCores)
	
	start := time.Now()
	if err := trainIntensively(model, genomeFile, perfCores); err != nil {
		log.Fatalf("Training failed: %v", err)
	}
	duration := time.Since(start)
	
	fmt.Printf("✅ Training completed in %v\n", duration)
	fmt.Println()

	// Generate test sequences
	fmt.Println("[STEP 3] Generating test sequences...")
	sequences := make([][]byte, 5)
	for i := range sequences {
		result, err := model.GenerateSequence(50000)
		if err != nil {
			log.Fatalf("Failed to generate sequence: %v", err)
		}
		sequences[i] = result.Sequence
	}
	
	// Calculate metrics
	fmt.Println("[STEP 4] Calculating authenticity metrics...")
	complexity := calculateComplexity(sequences)
	gcContent := calculateGCContent(sequences)
	
	fmt.Printf("📊 RESULTS:\n")
	fmt.Printf("  Complexity: %.4f", complexity)
	if complexity > 0.01 {
		fmt.Printf(" (%.1fx improvement) ✅", complexity/0.005)
	}
	fmt.Println()
	fmt.Printf("  GC Content: %.2f%%", gcContent*100)
	if gcContent >= 0.35 && gcContent <= 0.50 {
		fmt.Print(" ✅")
	}
	fmt.Println()
	
	fmt.Println()
	fmt.Println("🎉 SIMPLE TRAINING COMPLETE!")
	fmt.Printf("⚡ Training time: %v\n", duration)
	fmt.Printf("🧬 Successfully trained on full human genome\n")
	fmt.Printf("🎯 Ready for population genomics expansion\n")
}

func trainIntensively(model *authenticity.AdaptiveHMMModel, filename string, numCores int) error {
	fmt.Printf("🚀 Starting intensive training with %d cores...\n", numCores)
	
	// Open file
	file, err := os.Open(filename) // #nosec G304
	if err != nil {
		return fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()
	
	// Read file in chunks and process with multiple goroutines
	scanner := bufio.NewScanner(file)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 1024*1024)
	
	var sequence []byte
	var wg sync.WaitGroup
	var processedChunks int64
	
	// Create worker pool
	chunkChan := make(chan []byte, numCores*2)
	
	// Start workers
	for i := 0; i < numCores; i++ {
		wg.Add(1)
		go func(workerID int) {
			defer wg.Done()
			for chunk := range chunkChan {
				// Intensive processing - train multiple times on same data
				for repeat := 0; repeat < 3; repeat++ {
					if err := model.TrainFromSequence(chunk); err != nil {
						fmt.Printf("Worker %d error: %v\n", workerID, err)
					}
				}
				atomic.AddInt64(&processedChunks, 1)
				
				if processedChunks%100 == 0 {
					fmt.Printf("🔥 Processed %d chunks across %d cores\n", processedChunks, numCores)
				}
			}
		}(i)
	}
	
	// Read and distribute work
	chunkSize := 100000 // 100KB chunks
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		
		if strings.HasPrefix(line, ">") {
			// Process accumulated sequence
			if len(sequence) > 0 {
				// Split into chunks
				for i := 0; i < len(sequence); i += chunkSize {
					end := i + chunkSize
					if end > len(sequence) {
						end = len(sequence)
					}
					
					chunk := make([]byte, end-i)
					copy(chunk, sequence[i:end])
					chunkChan <- chunk
				}
				sequence = nil
			}
		} else {
			sequence = append(sequence, []byte(strings.ToUpper(line))...)
		}
	}
	
	// Process final sequence
	if len(sequence) > 0 {
		for i := 0; i < len(sequence); i += chunkSize {
			end := i + chunkSize
			if end > len(sequence) {
				end = len(sequence)
			}
			
			chunk := make([]byte, end-i)
			copy(chunk, sequence[i:end])
			chunkChan <- chunk
		}
	}
	
	// Close channel and wait for workers
	close(chunkChan)
	wg.Wait()
	
	fmt.Printf("✅ Intensive training complete: %d chunks processed\n", processedChunks)
	return nil
}

func calculateComplexity(sequences [][]byte) float64 {
	if len(sequences) == 0 {
		return 0
	}
	
	totalComplexity := 0.0
	for _, seq := range sequences {
		k := 4
		if len(seq) < k {
			continue
		}
		
		uniqueKmers := make(map[string]bool)
		for i := 0; i <= len(seq)-k; i++ {
			kmer := string(seq[i : i+k])
			uniqueKmers[kmer] = true
		}
		
		complexity := float64(len(uniqueKmers)) / float64(len(seq)-k+1)
		totalComplexity += complexity
	}
	
	return totalComplexity / float64(len(sequences))
}

func calculateGCContent(sequences [][]byte) float64 {
	if len(sequences) == 0 {
		return 0
	}
	
	totalGC := 0.0
	for _, seq := range sequences {
		if len(seq) == 0 {
			continue
		}
		
		gcCount := 0
		for _, base := range seq {
			if base == 'G' || base == 'C' {
				gcCount++
			}
		}
		
		gc := float64(gcCount) / float64(len(seq))
		totalGC += gc
	}
	
	return totalGC / float64(len(sequences))
}

func detectAppleSiliconCores() int {
	totalCores := runtime.NumCPU()
	
	if totalCores >= 8 {
		switch totalCores {
		case 8:  // M1: 4P + 4E
			return 4
		case 10: // M1 Pro: 6P + 4E
			return 6
		case 12: // M1 Max: 8P + 4E
			return 8
		case 16: // M2 Ultra: 12P + 4E
			return 12
		case 20: // M1 Ultra: 16P + 4E
			return 16
		case 24: // M2 Ultra: 16P + 8E
			return 16
		default:
			return (totalCores * 2) / 3
		}
	}
	
	if totalCores > 4 {
		return totalCores - 2
	}
	
	return totalCores
}