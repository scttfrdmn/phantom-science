// Command validate_vcf_sequences runs "go run ./cmd/experiments/sequence_validator"
// on sequences/validation (run it from the repository root), then tries to
// run a VCF validator that does not exist, and prints a fixed summary.
//
// It is a one-off experiment, not part of the supported phantom-science CLI.
// The VCF step always fails (vcf_validator.go was never written) and the
// comparison summary is hard-coded text, not a result of the run.
package main

import (
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"
)

func main() {
	log.SetFlags(log.LstdFlags)

	fmt.Println("=== VCF-TRAINED SEQUENCE VALIDATION ===")
	fmt.Println("Validating VCF-trained sequences vs synthetic-trained sequences")
	fmt.Println()

	// Validate synthetic sequences (original)
	fmt.Println("🧬 VALIDATING SYNTHETIC-TRAINED SEQUENCES:")
	fmt.Println("==========================================")

	cmd := exec.Command("go", "run", "./cmd/experiments/sequence_validator")
	cmd.Env = append(os.Environ(), "INPUT_DIR=sequences/validation", "OUTPUT_FILE=validation_synthetic.json")
	output, err := cmd.CombinedOutput()
	if err != nil {
		log.Printf("Warning: Failed to validate synthetic sequences: %v", err)
	} else {
		fmt.Print(string(output))
	}

	fmt.Println("\n" + strings.Repeat("=", 80))

	// Validate VCF-trained sequences
	fmt.Println("\n🧬 VALIDATING VCF-TRAINED SEQUENCES:")
	fmt.Println("====================================")

	// Create temporary validator for VCF sequences
	err = createVCFValidator()
	if err != nil {
		log.Printf("Failed to create VCF validator: %v", err)
		return
	}

	// NOTE: vcf_validator.go has never existed in this repository, so this step
	// always fails and is reported as a warning below.
	cmd = exec.Command("go", "run", "cmd/experiments/validate_vcf_sequences/vcf_validator.go")
	output, err = cmd.CombinedOutput()
	if err != nil {
		log.Printf("Warning: Failed to validate VCF sequences: %v", err)
	} else {
		fmt.Print(string(output))
	}

	// Compare results
	fmt.Println("\n" + strings.Repeat("=", 80))
	fmt.Println("\n📊 COMPARISON SUMMARY:")
	fmt.Println("======================")

	fmt.Println("Note: the text below is fixed, written after an earlier run; it is not")
	fmt.Println("derived from the output above.")
	fmt.Println()
	fmt.Println("⚠️  EARLIER FINDING: Both synthetic and VCF-trained models")
	fmt.Println("    are generating homopolymer sequences for certain populations")
	fmt.Println()
	fmt.Println("🔍 ROOT CAUSE ANALYSIS:")
	fmt.Println("   • Issue appears to be in the HMM model architecture itself")
	fmt.Println("   • Some populations (CEU, LWK, GWD, CHS, MXL) consistently fail")
	fmt.Println("   • Successful populations (YRI, CHB, JPT, PUR, TSI) show good diversity")
	fmt.Println()
	fmt.Println("🎯 NEXT STEPS:")
	fmt.Println("   1. Focus validation on the 5 working population models")
	fmt.Println("   2. Debug HMM model for homopolymer generation issue")
	fmt.Println("   3. Proceed with expert validation using successful models")
	fmt.Println("   4. Document model limitations in validation report")
}

func createVCFValidator() error {
	// This is a simplified approach - in practice, we'd modify the existing validator
	// to handle both input directories
	return nil
}
