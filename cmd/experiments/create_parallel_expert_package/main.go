// Command create_parallel_expert_package is the 11-population variant of
// create_expert_package: it copies the FASTA files under the hard-coded
// directory sequences/parallel_vcf_trained (the output of
// parallel_vcf_training) into expert_validation_package_parallel_vcf/ as
// header-stripped blind copies and labelled reference copies, and writes a
// Markdown evaluation protocol and README for reviewers. It also copies
// PARALLEL_VCF_TRAINING_BREAKTHROUGH.md from the working directory and exits
// with an error if that file is missing; the repository does not contain it.
// As in create_expert_package, the ID-to-file key is written into the
// package itself.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
package main

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"
	"time"
)

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== PARALLEL VCF EXPERT VALIDATION PACKAGE CREATION ===")
	fmt.Println("Creating anonymized parallel VCF-trained sequences for expert validation")
	fmt.Println("Using all 11 populations (models trained on 1000 Genomes VCF allele strings)")
	fmt.Println()

	// Source: parallel VCF-trained sequences from all populations
	sourceDir := "sequences/parallel_vcf_trained"
	outputDir := "expert_validation_package_parallel_vcf"

	// All 11 populations successfully trained
	allPopulations := []string{"YRI", "LWK", "GWD", "CEU", "GBR", "TSI", "CHB", "JPT", "CHS", "MXL", "PUR"}

	fmt.Printf("📁 Source directory: %s\n", sourceDir)
	fmt.Printf("📦 Output package: %s\n", outputDir)
	fmt.Printf("🧬 All populations: %v\n", allPopulations)
	fmt.Printf("🎯 Expected: 15 sequences per population\n")
	fmt.Println()

	// Create output directory
	err := os.MkdirAll(outputDir, 0750)
	if err != nil {
		log.Fatalf("Failed to create output directory: %v", err)
	}

	// Create package structure
	err = createPackageStructure(outputDir)
	if err != nil {
		log.Fatalf("Failed to create package structure: %v", err)
	}

	// Create anonymized sequences for blind testing
	totalSequences, err := createParallelAnonymizedSequences(sourceDir, outputDir, allPopulations)
	if err != nil {
		log.Fatalf("Failed to create anonymized sequences: %v", err)
	}

	// Create expert evaluation materials
	err = createParallelExpertMaterials(outputDir, totalSequences)
	if err != nil {
		log.Fatalf("Failed to create expert materials: %v", err)
	}

	fmt.Printf("\n✅ Parallel VCF expert validation package created successfully\n")
	fmt.Printf("📦 Package location: %s\n", outputDir)
	fmt.Printf("🧬 Total sequences: %d from %d populations\n", totalSequences, len(allPopulations))
	fmt.Printf("🔬 Package ready for review\n")
}

func createPackageStructure(outputDir string) error {
	subdirs := []string{
		"blind_test_sequences",
		"reference_sequences",
		"evaluation_materials",
		"results",
		"training_documentation",
	}

	for _, subdir := range subdirs {
		err := os.MkdirAll(filepath.Join(outputDir, subdir), 0750)
		if err != nil {
			return err
		}
	}

	fmt.Println("📁 Package structure created:")
	for _, subdir := range subdirs {
		fmt.Printf("   • %s/\n", subdir)
	}

	return nil
}

func createParallelAnonymizedSequences(sourceDir, outputDir string, populations []string) (int, error) {
	fmt.Println("\n🎭 Creating anonymized parallel VCF sequences for blind testing:")

	blindDir := filepath.Join(outputDir, "blind_test_sequences")
	refDir := filepath.Join(outputDir, "reference_sequences")

	// Create mapping file for de-anonymization
	mappingFile, err := os.Create(filepath.Join(outputDir, "sequence_mapping_CONFIDENTIAL.txt")) // #nosec G304
	if err != nil {
		return 0, err
	}
	defer mappingFile.Close()

	// Write mapping header
	mappingFile.WriteString("# CONFIDENTIAL: Parallel VCF Sequence Anonymization Mapping\n")
	mappingFile.WriteString("# Generated: " + time.Now().Format("2006-01-02 15:04:05") + "\n")
	mappingFile.WriteString("# Training: parallel_vcf_training (1000 Genomes VCF allele strings)\n")
	mappingFile.WriteString("# Data: 1000 Genomes Project VCF files (11 populations × 22 chromosomes)\n")
	mappingFile.WriteString("# Anonymous_ID -> Original_Population_File\n\n")

	totalSequences := 0
	anonymousID := 1

	// Process each population
	for _, pop := range populations {
		popDir := filepath.Join(sourceDir, pop)

		// Check if population directory exists
		if _, err := os.Stat(popDir); os.IsNotExist(err) {
			fmt.Printf("   ⚠️  %s: Directory not found, skipping\n", pop)
			continue
		}

		// Get all FASTA files for this population
		pattern := filepath.Join(popDir, "*.fasta")
		files, err := filepath.Glob(pattern)
		if err != nil {
			return 0, err
		}

		popSequences := 0

		// Create anonymized copies
		for _, file := range files {
			// Generate anonymous filename
			anonymousFilename := fmt.Sprintf("PARALLEL_VCF_SEQ_%04d.fasta", anonymousID)

			// Copy to blind test directory with anonymous name
			err = copyAndAnonymizeParallelFile(file, filepath.Join(blindDir, anonymousFilename))
			if err != nil {
				log.Printf("Warning: Failed to anonymize %s: %v", filepath.Base(file), err)
				continue
			}

			// Copy to reference directory with population prefix
			err = copyFile(file, filepath.Join(refDir, pop+"_"+filepath.Base(file)))
			if err != nil {
				log.Printf("Warning: Failed to copy reference %s: %v", filepath.Base(file), err)
			}

			// Record mapping
			mappingFile.WriteString(fmt.Sprintf("PARALLEL_VCF_SEQ_%04d.fasta -> %s/%s\n",
				anonymousID, pop, filepath.Base(file)))

			anonymousID++
			popSequences++
			totalSequences++
		}

		fmt.Printf("   ✅ %s: %d parallel VCF sequences anonymized\n", pop, popSequences)
	}

	return totalSequences, nil
}

func copyAndAnonymizeParallelFile(srcPath, dstPath string) error {
	// Read original file
	content, err := os.ReadFile(srcPath) // #nosec G304
	if err != nil {
		return err
	}

	// Anonymize header (remove population-specific information)
	anonymizedContent := anonymizeParallelSequence(content)

	// Write anonymized file
	return os.WriteFile(dstPath, anonymizedContent, 0600)
}

func copyFile(srcPath, dstPath string) error {
	content, err := os.ReadFile(srcPath) // #nosec G304
	if err != nil {
		return err
	}
	return os.WriteFile(dstPath, content, 0600)
}

func anonymizeParallelSequence(content []byte) []byte {
	contentStr := string(content)
	lines := strings.Split(contentStr, "\n")

	if len(lines) == 0 {
		return content
	}

	// Anonymize the header line (first line starting with >)
	if strings.HasPrefix(lines[0], ">") {
		// Extract just length and basic quality metrics, remove population info
		length := extractMetadata(lines[0], "length")
		complexity := extractMetadata(lines[0], "complexity")
		gcContent := extractMetadata(lines[0], "gc_content")

		// Create anonymous header with key metrics
		lines[0] = fmt.Sprintf(">ANONYMOUS_PARALLEL_VCF_SEQUENCE length=%s complexity=%s gc_content=%s vcf_trained=true",
			length, complexity, gcContent)
	}

	return []byte(strings.Join(lines, "\n"))
}

func extractMetadata(header, key string) string {
	// Extract metadata from header
	fields := strings.Fields(header)
	for _, field := range fields {
		if strings.HasPrefix(field, key+"=") {
			return strings.TrimPrefix(field, key+"=")
		}
	}
	return "unknown"
}

func createParallelExpertMaterials(outputDir string, totalSequences int) error {
	fmt.Println("\n📋 Creating parallel VCF expert evaluation materials:")

	materialsDir := filepath.Join(outputDir, "evaluation_materials")

	// Create evaluation protocol
	err := createParallelEvaluationProtocol(materialsDir, totalSequences)
	if err != nil {
		return err
	}

	// Create expert README
	err = createParallelExpertREADME(materialsDir)
	if err != nil {
		return err
	}

	// Create training documentation
	err = createTrainingDocumentation(outputDir)
	if err != nil {
		return err
	}

	fmt.Println("   ✅ Parallel VCF evaluation protocol created")
	fmt.Println("   ✅ Expert README created")
	fmt.Println("   ✅ Training documentation included")

	return nil
}

func createParallelEvaluationProtocol(materialsDir string, totalSequences int) error {
	protocolPath := filepath.Join(materialsDir, "PARALLEL_VCF_EVALUATION_PROTOCOL.md")

	protocol := fmt.Sprintf(`# PARALLEL VCF GENOMIC SEQUENCE AUTHENTICITY EVALUATION PROTOCOL

## Overview

**Objective**: Assess the biological authenticity of parallel VCF-trained genomic sequences

**Method**: k=8 adaptive HMMs, one per population, trained by parallel_vcf_training on the REF and ALT allele strings of 1000 Genomes VCF records concatenated together (not on reference sequence)

**Sequences**: %d anonymized sequences from 11 global populations

**Status**: Experimental; the models have not been validated against real population data

---

## Training Background

### **Data Source**
- **1000 Genomes Project** VCF files, 22 autosomes per population
- Training input is the concatenated REF/ALT allele strings of the VCF records, so it has no reference-sequence context

### **Populations Included**
- **Africa**: YRI (Yoruba), LWK (Luhya), GWD (Gambian)
- **Europe**: CEU (Utah Europeans), GBR (British), TSI (Tuscans)  
- **East Asia**: CHB (Chinese), JPT (Japanese), CHS (Southern Chinese)
- **Americas**: MXL (Mexican), PUR (Puerto Rican)

### **Limitations**
- Sequence statistics for this package are not computed by this program; see the parallel_vcf_training output
- The population models have not been compared with real sequence from the same populations

---

## Evaluation Criteria

### 1. Biological Authenticity Assessment (Scale: 1-10)

**Score 9-10**: Indistinguishable from real genomic sequences
- Natural k-mer distribution patterns matching population genomics
- Appropriate GC content variation (47-52%% range)
- Realistic homopolymer runs (<=20 bases typical)
- Population-specific genomic signatures present

**Score 7-8**: High authenticity with subtle computational artifacts
- Mostly natural sequence characteristics
- Minor deviations in expected population patterns
- Generally consistent with real genomic data

**Score 5-6**: Moderate authenticity with noticeable issues
- Some unnatural sequence patterns detected
- Population signatures unclear or inconsistent
- Mixed authentic and artificial characteristics

**Score 3-4**: Low authenticity with significant problems
- Clear computational artifacts present
- Unnatural k-mer distributions
- Population characteristics not credible

**Score 1-2**: Poor authenticity, obviously artificial
- Extreme computational patterns (e.g., homopolymers)
- No biological credibility
- Unsuitable for any genomic applications

### 2. Population Origin Identification (Optional)

Can you identify population ancestry from sequence patterns?
- **African populations**: Higher genetic diversity, specific SNP patterns
- **European populations**: Lower diversity, founder effect signatures  
- **East Asian populations**: Distinct allele frequency patterns
- **American populations**: Admixed genomic characteristics
- **Cannot determine**: Population origin unclear from sequence

### 3. Comparative Assessment

How do these sequences compare to:
- Real human genomic sequences from public databases
- Other computational sequence generators (if familiar)
- Expected population genomic diversity patterns
- Clinical/research quality genomic data

---

## Evaluation Process

### Phase 1: Blind Assessment (Primary)
1. Review anonymized sequences in blind_test_sequences/
2. Score each sequence for biological authenticity (1-10)
3. Attempt population origin identification (optional)
4. Note any concerning patterns or obvious artifacts
5. **Do not** consult reference sequences during this phase

### Phase 2: Reference Comparison (Secondary)  
1. Review reference sequences with population labels
2. Assess population-specific pattern consistency
3. Validate or revise initial assessments
4. Compare your population identifications with actual labels

### Phase 3: Overall Assessment
1. Provide summary authenticity assessment
2. Discuss whether the sequences are fit for testing and research use
3. Compare with other genomic sequence sources
4. Recommend improvements or validation concerns

---

## Key Evaluation Questions

1. **Authenticity**: Do these sequences appear biologically realistic?
2. **Population Specificity**: Can you detect population-specific genomic patterns?
3. **Quality**: Are these suitable for genomic research applications?
4. **Artifacts**: Do you notice any computational generation artifacts?
5. **Utility**: How would these compare to real genomic sequences for research?

---

## Expected Timeline

- **Blind Assessment**: 2-3 hours for %d sequences
- **Reference Comparison**: 1-2 hours for validation
- **Report Preparation**: 1 hour for detailed feedback
- **Discussion Session**: 1 hour with research team

---

## Technical Context

The models were trained per population from 1000 Genomes VCF allele strings, processed in parallel across 11 populations. Your assessment will help determine how far the generated sequences are from real human sequence and whether population differences are detectable.

---

**Confidentiality**: All sequences and evaluation materials are confidential research data.
`, totalSequences, totalSequences)

	return os.WriteFile(protocolPath, []byte(protocol), 0600)
}

func createParallelExpertREADME(materialsDir string) error {
	readmePath := filepath.Join(materialsDir, "README_PARALLEL_VCF_EXPERT_VALIDATION.md")

	readme := `# PARALLEL VCF EXPERT VALIDATION PACKAGE - README

## Welcome, Genomics Expert!

Thank you for reviewing these computationally generated genomic sequences.

---

## Background

The sequences were generated by k=8 adaptive HMMs, one per population, for 11 1000 Genomes populations. Each model was trained on the REF and ALT allele strings of that population's VCF records concatenated together, not on reference sequence. The approach is experimental and has not been validated against real sequence.

---

## Package Contents

### 📁 **blind_test_sequences/**
Anonymized sequences for blind evaluation (primary assessment)
- File format: PARALLEL_VCF_SEQ_0001.fasta, etc.
- Headers anonymized to prevent bias
- Drawn from all 11 populations

### 📁 **reference_sequences/**
Same sequences with population labels for comparison
- File format: POPULATION_originalname.fasta
- Use after blind assessment for validation
- Enables population-specific pattern assessment

### 📁 **evaluation_materials/**
- **PARALLEL_VCF_EVALUATION_PROTOCOL.md**: Detailed evaluation guidelines
- **README_PARALLEL_VCF_EXPERT_VALIDATION.md**: This file

### 📁 **training_documentation/**
- **PARALLEL_VCF_TRAINING_BREAKTHROUGH.md**: training notes copied from the working directory (not reviewed for accuracy)

### 📁 **results/**
Empty directory for your completed evaluation forms

---

## Evaluation Process

### Step 1: Review Training Documentation
Read how the models were trained and what data they saw.

### Step 2: Blind Assessment (Critical)
1. Examine sequences in blind_test_sequences/
2. Score biological authenticity (1-10 scale)
3. Attempt population identification (optional)
4. **Do not** look at reference sequences yet

### Step 3: Reference Validation
1. Compare your assessments with population labels
2. Assess population-specific genomic patterns
3. Validate or revise initial scoring

### Step 4: Comprehensive Report
1. Provide overall authenticity assessment
2. Discuss whether the sequences are fit for testing and research use
3. Compare with real genomic sequences
4. Schedule discussion with research team

---

## Critical Evaluation Focus

### **Key Questions**
1. Do these sequences appear biologically authentic?
2. Can you detect population-specific genomic signatures?
3. Are these suitable for testing and research use?
4. How do they compare to real genomic sequences?
5. Do you observe any computational artifacts?

### **Population Genomics Context**
- **African populations**: Highest genetic diversity, ancient lineages
- **European populations**: Founder effects, reduced diversity
- **East Asian populations**: Distinct allele frequency patterns
- **American populations**: Recent admixture, mixed ancestry

---

## Expected Outcomes

### **Quality Benchmarks**
- Biological authenticity scores >7/10
- Detectable population-specific signatures
- Suitable for genomic research applications
- Comparable quality to real genomic sequences

---

## Important Notes

⚠️  **Confidentiality**: All materials are confidential research data

🎯 **Objectivity**: Please provide honest, critical assessment

🔬 **Expertise**: Your genomics knowledge is essential for validation

📊 **Impact**: Results will guide further development

---

**Thank you for your help.**
`

	return os.WriteFile(readmePath, []byte(readme), 0600)
}

func createTrainingDocumentation(outputDir string) error {
	trainingDir := filepath.Join(outputDir, "training_documentation")

	// Copy the training notes
	srcPath := "PARALLEL_VCF_TRAINING_BREAKTHROUGH.md"
	dstPath := filepath.Join(trainingDir, "PARALLEL_VCF_TRAINING_BREAKTHROUGH.md")

	content, err := os.ReadFile(srcPath)
	if err != nil {
		return fmt.Errorf("failed to read training documentation: %v", err)
	}

	return os.WriteFile(dstPath, content, 0600)
}
