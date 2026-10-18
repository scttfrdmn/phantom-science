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
	
	fmt.Println("=== EXPERT VALIDATION PACKAGE CREATION ===")
	fmt.Println("Creating anonymized sequences for expert validation")
	fmt.Println()
	
	// Source: validated sequences from working populations
	sourceDir := "sequences/validation"
	outputDir := "expert_validation_package"
	
	// Working populations (from validation analysis)
	workingPopulations := []string{"YRI", "CHB", "JPT", "PUR", "TSI"}
	
	fmt.Printf("📁 Source directory: %s\n", sourceDir)
	fmt.Printf("📦 Output package: %s\n", outputDir)
	fmt.Printf("🧬 Working populations: %v\n", workingPopulations)
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
	totalSequences, err := createAnonymizedSequences(sourceDir, outputDir, workingPopulations)
	if err != nil {
		log.Fatalf("Failed to create anonymized sequences: %v", err)
	}
	
	// Create expert evaluation materials
	err = createExpertMaterials(outputDir, totalSequences)
	if err != nil {
		log.Fatalf("Failed to create expert materials: %v", err)
	}
	
	fmt.Printf("\n✅ Expert validation package created successfully\n")
	fmt.Printf("📦 Package location: %s\n", outputDir)
	fmt.Printf("🧬 Total sequences: %d from %d populations\n", totalSequences, len(workingPopulations))
	fmt.Printf("🔬 Ready for expert validation\n")
}

func createPackageStructure(outputDir string) error {
	subdirs := []string{
		"blind_test_sequences",
		"reference_sequences", 
		"evaluation_materials",
		"results",
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

func createAnonymizedSequences(sourceDir, outputDir string, populations []string) (int, error) {
	fmt.Println("\n🎭 Creating anonymized sequences for blind testing:")
	
	blindDir := filepath.Join(outputDir, "blind_test_sequences")
	refDir := filepath.Join(outputDir, "reference_sequences")
	
	// Create mapping file for de-anonymization
	mappingFile, err := os.Create(filepath.Join(outputDir, "sequence_mapping_CONFIDENTIAL.txt")) // #nosec G304
	if err != nil {
		return 0, err
	}
	defer mappingFile.Close()
	
	// Write mapping header
	mappingFile.WriteString("# CONFIDENTIAL: Sequence Anonymization Mapping\n")
	mappingFile.WriteString("# Generated: " + time.Now().Format("2006-01-02 15:04:05") + "\n")
	mappingFile.WriteString("# Anonymous_ID -> Original_Population_File\n\n")
	
	totalSequences := 0
	anonymousID := 1
	
	
	// Process each working population
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
			anonymousFilename := fmt.Sprintf("SEQ_%04d.fasta", anonymousID)
			
			// Copy to blind test directory with anonymous name
			err = copyAndAnonymizeFile(file, filepath.Join(blindDir, anonymousFilename))
			if err != nil {
				log.Printf("Warning: Failed to anonymize %s: %v", filepath.Base(file), err)
				continue
			}
			
			// Copy to reference directory with original name for expert reference
			err = copyFile(file, filepath.Join(refDir, pop+"_"+filepath.Base(file)))
			if err != nil {
				log.Printf("Warning: Failed to copy reference %s: %v", filepath.Base(file), err)
			}
			
			// Record mapping
			mappingFile.WriteString(fmt.Sprintf("SEQ_%04d.fasta -> %s/%s\n", 
				anonymousID, pop, filepath.Base(file)))
			
			anonymousID++
			popSequences++
			totalSequences++
		}
		
		fmt.Printf("   ✅ %s: %d sequences anonymized\n", pop, popSequences)
	}
	
	return totalSequences, nil
}

func copyAndAnonymizeFile(srcPath, dstPath string) error {
	// Read original file
	content, err := os.ReadFile(srcPath) // #nosec G304
	if err != nil {
		return err
	}
	
	// Anonymize header (remove population-specific information)
	anonymizedContent := anonymizeSequence(content)
	
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

func anonymizeSequence(content []byte) []byte {
	contentStr := string(content)
	lines := strings.Split(contentStr, "\n")
	
	if len(lines) == 0 {
		return content
	}
	
	// Anonymize the header line (first line starting with >)
	if strings.HasPrefix(lines[0], ">") {
		// Extract just the sequence ID, remove population/region info
		parts := strings.Fields(lines[0])
		if len(parts) > 0 {
			// Create anonymous header
			lines[0] = ">ANONYMOUS_SEQUENCE length=" + extractLength(lines[0])
		}
	}
	
	return []byte(strings.Join(lines, "\n"))
}

func extractLength(header string) string {
	// Extract length from header metadata
	fields := strings.Fields(header)
	for _, field := range fields {
		if strings.HasPrefix(field, "length=") {
			return strings.TrimPrefix(field, "length=")
		}
	}
	return "unknown"
}

func createExpertMaterials(outputDir string, totalSequences int) error {
	fmt.Println("\n📋 Creating expert evaluation materials:")
	
	materialsDir := filepath.Join(outputDir, "evaluation_materials")
	
	// Create evaluation protocol
	err := createEvaluationProtocol(materialsDir, totalSequences)
	if err != nil {
		return err
	}
	
	// Create expert README
	err = createExpertREADME(materialsDir)
	if err != nil {
		return err
	}
	
	fmt.Println("   ✅ Evaluation protocol created")
	fmt.Println("   ✅ Expert README created")
	
	return nil
}

func createEvaluationProtocol(materialsDir string, totalSequences int) error {
	protocolPath := filepath.Join(materialsDir, "EVALUATION_PROTOCOL.md")
	
	protocol := fmt.Sprintf(`# GENOMIC SEQUENCE AUTHENTICITY EVALUATION PROTOCOL

## Overview

**Objective**: Assess the biological authenticity of computationally generated genomic sequences

**Method**: Population-specific k-mer Hidden Markov Models (k=8)

**Sequences**: %d anonymized sequences from 5 global populations

**Evaluation**: Blind assessment by genomics specialists

---

## Evaluation Criteria

### 1. Biological Authenticity (Scale: 1-10)

**Score 9-10**: Indistinguishable from real genomic sequences
- Natural k-mer distribution patterns
- Appropriate GC content variation
- Realistic homopolymer runs
- Biologically plausible sequence characteristics

**Score 7-8**: High authenticity with minor artificial patterns
- Mostly natural sequence characteristics
- Slight deviations in expected patterns
- Generally biologically credible

**Score 5-6**: Moderate authenticity with noticeable artificial elements
- Some unnatural sequence patterns
- Detectable computational artifacts
- Mixed authentic and artificial characteristics

**Score 3-4**: Low authenticity with significant artificial patterns
- Clear computational artifacts
- Unnatural sequence distributions
- Limited biological credibility

**Score 1-2**: Poor authenticity, clearly artificial
- Obvious computational generation patterns
- Biologically implausible characteristics
- Not suitable for biological applications

### 2. Population Origin Identification (Optional)

If sequence patterns suggest specific population characteristics:
- African populations (high genetic diversity)
- European populations (reduced diversity, founder effects)
- East Asian populations (distinct genetic signatures)
- American populations (admixed characteristics)
- Unable to determine population origin

### 3. Technical Quality Assessment

**Sequence Composition**:
- GC content appropriateness
- K-mer complexity and diversity
- Homopolymer run lengths
- Overall sequence quality

**Comparative Assessment**:
- How do these sequences compare to:
  - Real human genomic sequences
  - Other computational sequence generators
  - Expected population genomic patterns

---

## Evaluation Process

### Phase 1: Blind Assessment
1. Review anonymized sequences (blind_test_sequences/)
2. Complete evaluation for each sequence
3. Provide overall authenticity scores
4. Note any concerning patterns or artifacts

### Phase 2: Reference Comparison
1. Review reference sequences with population labels
2. Assess population-specific patterns
3. Validate or revise initial assessments

### Phase 3: Discussion and Feedback
1. Discuss findings with research team
2. Provide recommendations for improvement
3. Assess potential clinical/research applications

---

## Expected Timeline

- **Initial Review**: 2-3 hours for sequence assessment
- **Detailed Analysis**: 1-2 hours for comparative analysis
- **Discussion Session**: 1 hour with research team

---

**Confidentiality**: All sequences and evaluation materials are confidential research data.
`, totalSequences)

	return os.WriteFile(protocolPath, []byte(protocol), 0600)
}

func createExpertREADME(materialsDir string) error {
	readmePath := filepath.Join(materialsDir, "README_FOR_EXPERTS.md")
	
	readme := `# EXPERT VALIDATION PACKAGE - README

## Welcome, Genomics Specialist!

Thank you for participating in this important validation study of computationally generated genomic sequences.

---

## Package Contents

### 📁 **blind_test_sequences/**
Anonymized sequences for blind evaluation. These are the primary sequences for your assessment.
- File format: SEQ_0001.fasta, SEQ_0002.fasta, etc.
- Headers anonymized to prevent bias
- Generated from population-specific k=8 HMM models

### 📁 **reference_sequences/**
Same sequences with population labels for reference comparison.
- File format: POPULATION_originalname.fasta
- Use after blind assessment for validation
- Helps assess population-specific patterns

### 📁 **evaluation_materials/**
- **EVALUATION_PROTOCOL.md**: Detailed evaluation guidelines
- **README_FOR_EXPERTS.md**: This file

### 📁 **results/**
Empty directory for your completed evaluation forms.

---

## Evaluation Process

### Step 1: Read Protocol
Review EVALUATION_PROTOCOL.md carefully for:
- Scoring criteria (1-10 scale)
- Assessment guidelines
- Technical quality metrics

### Step 2: Blind Assessment
1. Examine sequences in blind_test_sequences/
2. Assess biological authenticity and technical quality
3. **Do not** look at reference sequences yet

### Step 3: Reference Comparison
1. Review your assessments
2. Compare with population-labeled references
3. Note any patterns or revisions

### Step 4: Submit Results
1. Schedule discussion session with research team

---

## Research Context

**Project**: Population-specific genomic sequence generation using adaptive k-mer HMM models

**Innovation**: First approach to generate population-authentic genomic sequences for precision medicine applications

**Your Role**: Critical validation of biological authenticity and clinical suitability

---

## Technical Background

**Method**: k=8 Hidden Markov Models trained on population-specific genomic data

**Populations**: 5 global populations (anonymized during blind testing)

**Metrics**: K-mer complexity, GC content, homopolymer analysis, biological authenticity

**Goal**: Generate sequences indistinguishable from real human genomic sequences

---

## Important Notes

⚠️  **Confidentiality**: All materials are confidential research data

🎯 **Objectivity**: Please provide honest, critical assessment

🔬 **Expertise**: Your genomics knowledge is essential for validation

📊 **Impact**: Results will guide clinical applications and further development

---

**Thank you for your valuable contribution to advancing population genomics!**
`

	return os.WriteFile(readmePath, []byte(readme), 0600)
}