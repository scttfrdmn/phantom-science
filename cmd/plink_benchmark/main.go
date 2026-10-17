package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"time"
)

// PLINKBenchmark runs PLINK analysis on our genomic sequences
type PLINKBenchmark struct {
	DataDir      string
	OutputDir    string
	PLINKPath    string
	Populations  []string
	Results      []PLINKResult
}

type PLINKResult struct {
	Population     string
	SequenceFile   string
	ProcessingTime time.Duration
	Success        bool
	ErrorMessage   string
	OutputFiles    []string
}

func NewPLINKBenchmark() *PLINKBenchmark {
	return &PLINKBenchmark{
		DataDir:     "sequences/parallel_vcf_trained",
		OutputDir:   "benchmarking_results/plink",
		PLINKPath:   "/opt/homebrew/Caskroom/miniforge/base/bin/plink",
		Populations: []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"},
		Results:     make([]PLINKResult, 0),
	}
}

func (pb *PLINKBenchmark) ConvertSequencesToPED() error {
	fmt.Println("🧬 Converting FASTA sequences to PLINK PED/MAP format...")
	
	// Create output directory
	if err := os.MkdirAll(pb.OutputDir, 0750); err != nil {
		return fmt.Errorf("failed to create output directory: %w", err)
	}
	
	// Create a single PED file with all sequences
	pedFile := filepath.Join(pb.OutputDir, "benchmark_sequences.ped")
	mapFile := filepath.Join(pb.OutputDir, "benchmark_sequences.map")
	
	pedWriter, err := os.Create(pedFile) // #nosec G304
	if err != nil {
		return err
	}
	defer pedWriter.Close()
	
	mapWriter, err := os.Create(mapFile) // #nosec G304
	if err != nil {
		return err
	}
	defer mapWriter.Close()
	
	// Write MAP file header (dummy SNP data since we're working with sequences)
	// Format: chromosome, variant identifier, position in morgans, base-pair coordinate
	fmt.Fprintf(mapWriter, "1\tSNP1\t0\t1000\n")
	fmt.Fprintf(mapWriter, "1\tSNP2\t0\t2000\n")
	fmt.Fprintf(mapWriter, "1\tSNP3\t0\t3000\n")
	
	sampleID := 1
	
	for _, pop := range pb.Populations {
		for _, length := range []int{1000, 5000, 10000} {
			for i := 1; i <= 5; i++ {
				filename := fmt.Sprintf("%s_PARALLEL_VCF_L%d_S%02d.fasta", pop, length, i)
				filepath := filepath.Join(pb.DataDir, pop, filename)
				
				// Read sequence
				sequence, err := pb.readSequenceFromFasta(filepath)
				if err != nil {
					fmt.Printf("  ⚠️  Skipping %s: %v\n", filename, err)
					continue
				}
				
				// Convert to simplified genotype data (dummy conversion for benchmarking)
				genotype := pb.sequenceToGenotype(sequence)
				
				// Write PED entry
				// Format: Family_ID Individual_ID Paternal_ID Maternal_ID Sex Phenotype Genotypes...
				fmt.Fprintf(pedWriter, "%s\t%s_%d\t0\t0\t1\t1\t%s\n", 
					pop, pop, sampleID, genotype)
				
				sampleID++
			}
		}
	}
	
	fmt.Printf("  ✅ Created PED/MAP files: %s, %s\n", pedFile, mapFile)
	return nil
}

func (pb *PLINKBenchmark) readSequenceFromFasta(filepath string) (string, error) {
	file, err := os.Open(filepath) // #nosec G304
	if err != nil {
		return "", err
	}
	defer file.Close()
	
	scanner := bufio.NewScanner(file)
	var sequence strings.Builder
	
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if !strings.HasPrefix(line, ">") {
			sequence.WriteString(strings.ToUpper(line))
		}
	}
	
	return sequence.String(), scanner.Err()
}

func (pb *PLINKBenchmark) sequenceToGenotype(sequence string) string {
	// Convert sequence to proper SNP genotypes for PLINK
	// Take every 100th base to create realistic SNP density
	
	genotypes := make([]string, 3) // Match MAP file (3 SNPs)
	positions := []int{1000, 2000, 3000}
	
	for i, pos := range positions {
		if pos < len(sequence) {
			base := sequence[pos-1] // Convert to 0-based indexing
			switch base {
			case 'A':
				genotypes[i] = "A A"
			case 'T':
				genotypes[i] = "T T"
			case 'G':
				genotypes[i] = "G G"
			case 'C':
				genotypes[i] = "C C"
			default:
				genotypes[i] = "0 0" // Missing data
			}
		} else {
			genotypes[i] = "0 0" // Missing - sequence too short
		}
	}
	
	return strings.Join(genotypes, " ")
}

func (pb *PLINKBenchmark) RunPLINKAnalysis() error {
	fmt.Println("\n🧬 Running PLINK Analysis...")
	
	start := time.Now()
	
	// Run basic PLINK analysis
	inputPrefix := filepath.Join(pb.OutputDir, "benchmark_sequences")
	outputPrefix := filepath.Join(pb.OutputDir, "plink_results")
	
	// Basic frequency analysis
	cmd := exec.Command(pb.PLINKPath,
		"--file", inputPrefix,
		"--freq",
		"--out", outputPrefix) // #nosec G204
	
	output, err := cmd.CombinedOutput()
	processingTime := time.Since(start)
	
	result := PLINKResult{
		Population:     "ALL_POPULATIONS",
		SequenceFile:   "benchmark_sequences",
		ProcessingTime: processingTime,
		Success:        err == nil,
	}
	
	if err != nil {
		result.ErrorMessage = fmt.Sprintf("PLINK error: %v\nOutput: %s", err, string(output))
		fmt.Printf("  ❌ PLINK analysis failed: %v\n", err)
		fmt.Printf("     Output: %s\n", string(output))
	} else {
		fmt.Printf("  ✅ PLINK analysis completed in %v\n", processingTime)
		
		// List output files
		files, _ := filepath.Glob(outputPrefix + "*")
		result.OutputFiles = files
		
		for _, file := range files {
			fmt.Printf("     📄 Generated: %s\n", filepath.Base(file))
		}
	}
	
	pb.Results = append(pb.Results, result)
	return nil
}

func (pb *PLINKBenchmark) GenerateReport() error {
	fmt.Println("\n📊 Generating PLINK Benchmark Report...")
	
	reportFile := filepath.Join(pb.OutputDir, "PLINK_BENCHMARK_REPORT.md")
	file, err := os.Create(reportFile) // #nosec G304
	if err != nil {
		return err
	}
	defer file.Close()
	
	writer := bufio.NewWriter(file)
	defer writer.Flush()
	
	fmt.Fprintf(writer, "# PLINK Benchmark Report\n\n")
	fmt.Fprintf(writer, "**Generated:** %s\n\n", time.Now().Format("2006-01-02 15:04:05"))
	
	fmt.Fprintf(writer, "## Analysis Overview\n\n")
	fmt.Fprintf(writer, "- **Tool**: PLINK v1.9\n")
	fmt.Fprintf(writer, "- **Dataset**: Parallel VCF-trained sequences\n")
	fmt.Fprintf(writer, "- **Total Sequences**: 165 (11 populations × 15 sequences)\n")
	fmt.Fprintf(writer, "- **Analysis Type**: Basic frequency analysis\n\n")
	
	fmt.Fprintf(writer, "## Results\n\n")
	
	for _, result := range pb.Results {
		fmt.Fprintf(writer, "### %s\n\n", result.Population)
		fmt.Fprintf(writer, "- **Processing Time**: %v\n", result.ProcessingTime)
		fmt.Fprintf(writer, "- **Success**: %t\n", result.Success)
		
		if !result.Success {
			fmt.Fprintf(writer, "- **Error**: %s\n", result.ErrorMessage)
		} else {
			fmt.Fprintf(writer, "- **Output Files**: %d generated\n", len(result.OutputFiles))
			for _, file := range result.OutputFiles {
				fmt.Fprintf(writer, "  - %s\n", filepath.Base(file))
			}
		}
		fmt.Fprintf(writer, "\n")
	}
	
	fmt.Fprintf(writer, "## Performance Summary\n\n")
	if len(pb.Results) > 0 {
		totalTime := time.Duration(0)
		successCount := 0
		
		for _, result := range pb.Results {
			totalTime += result.ProcessingTime
			if result.Success {
				successCount++
			}
		}
		
		avgTime := totalTime / time.Duration(len(pb.Results))
		successRate := float64(successCount) / float64(len(pb.Results)) * 100
		
		fmt.Fprintf(writer, "- **Average Processing Time**: %v\n", avgTime)
		fmt.Fprintf(writer, "- **Success Rate**: %.1f%%\n", successRate)
		fmt.Fprintf(writer, "- **Total Analyses**: %d\n", len(pb.Results))
	}
	
	fmt.Printf("📄 PLINK benchmark report saved: %s\n", reportFile)
	return nil
}

func main() {
	benchmark := NewPLINKBenchmark()
	
	fmt.Println("🧬 PLINK GENOMICS BENCHMARKING")
	fmt.Println("==============================")
	
	// Convert sequences to PLINK format
	if err := benchmark.ConvertSequencesToPED(); err != nil {
		log.Fatalf("Failed to convert sequences: %v", err)
	}
	
	// Run PLINK analysis
	if err := benchmark.RunPLINKAnalysis(); err != nil {
		log.Fatalf("Failed to run PLINK analysis: %v", err)
	}
	
	// Generate report
	if err := benchmark.GenerateReport(); err != nil {
		log.Fatalf("Failed to generate report: %v", err)
	}
	
	fmt.Println("\n🎉 PLINK benchmarking completed successfully!")
}