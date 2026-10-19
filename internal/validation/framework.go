// Package validation runs simple checks on FASTQ data: record structure,
// quality-score range and encoding, read lengths, base composition, GC
// content, homopolymer reads, and a Hardy-Weinberg chi-squared test.
//
// The checks catch malformed or degenerate output; they do not show that
// data is biologically realistic. Each result's Confidence is a fixed value
// chosen per test (the HWE test reports its p-value), not a statistical
// confidence. The ComplianceLevel and RegulatoryReady fields of the results
// only echo configuration flags; no regulatory check is performed.
package validation

import (
	"fmt"
	"math"
	"slices"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/authenticity/statistics"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// GenotypeCount holds observed genotype counts for a single variant.
type GenotypeCount struct {
	NHomRef int // Homozygous reference (AA)
	NHet    int // Heterozygous (Aa)
	NHomAlt int // Homozygous alternate (aa)
}

// ValidationData carries input for the validation framework.
// FASTQ is used for format and sequence-based tests.
// GenotypeCounts, if non-empty, provide explicit genotype observations for HWE tests
// instead of the GC-proxy derivation.
type ValidationData struct {
	FASTQ          string
	GenotypeCounts []GenotypeCount
}

// NewValidationData creates ValidationData wrapping a FASTQ string.
func NewValidationData(fastq string) *ValidationData {
	return &ValidationData{FASTQ: fastq}
}

// Framework runs the tests named in a types.ValidationConfig.
type Framework struct {
	config      types.ValidationConfig
	readLengths ReadLengthExpectation
}

// Option configures a Framework.
type Option func(*Framework)

// NewFramework returns a framework running the tests listed in config.
// Without options, read_length_consistency requires every read to have the
// same length (ReadLengthsUniform).
func NewFramework(config types.ValidationConfig, opts ...Option) *Framework {
	f := &Framework{config: config}
	for _, opt := range opts {
		opt(f)
	}
	return f
}

// ReadLengthMode selects how the read_length_consistency test judges read
// lengths.
type ReadLengthMode int

const (
	// ReadLengthsUniform (the default) requires every read to have the same
	// length, and that length to be Length when Length is set. Use it for
	// fixed-length (Illumina) data.
	ReadLengthsUniform ReadLengthMode = iota
	// ReadLengthsAny accepts any lengths and reports their summary
	// (count, min, max, mean, N50) as an informational pass. Use it when the
	// expected distribution is unknown.
	ReadLengthsAny
	// ReadLengthsDistributed accepts variable lengths within [Min, Max] and,
	// when Mean and SD are set, requires the sample mean of n reads to lie
	// within 4·SD/√n of Mean.
	ReadLengthsDistributed
)

// ReadLengthExpectation describes the read lengths a FASTQ file should have.
type ReadLengthExpectation struct {
	Mode ReadLengthMode
	// Length is the required read length for ReadLengthsUniform; 0 accepts
	// any common length.
	Length int
	// Min and Max are inclusive bounds for ReadLengthsDistributed; 0 leaves
	// that side unchecked.
	Min, Max int
	// Mean and SD are the expected mean and per-read standard deviation for
	// ReadLengthsDistributed; the mean is checked only when both are > 0.
	Mean, SD float64
	// Name describes the distribution in messages (e.g. "lognormal").
	Name string
}

// WithReadLengths sets the expectation used by read_length_consistency.
func WithReadLengths(exp ReadLengthExpectation) Option {
	return func(f *Framework) { f.readLengths = exp }
}

// RunValidation runs every configured test against data. AcademicReady is
// set only when all tests ran and passed and academic validation
// (PeerReviewReady) was requested; RegulatoryReady is never set, because no
// regulatory validation procedure exists. See types.ValidationResults.
func (f *Framework) RunValidation(data interface{}) types.ValidationResults {
	// Check if validation is disabled
	if !f.config.Enable {
		return types.ValidationResults{
			OverallPass:     false,
			ComplianceLevel: "none",
			Tests:           []types.ValidationResult{},
			Summary:         "Validation is disabled",
			RegulatoryReady: false,
			AcademicReady:   false,
		}
	}

	// Readiness is decided after the tests run. RegulatoryReady stays false:
	// no FDA/EMA validation procedure is implemented, so the regulatory
	// config flags cannot support a regulatory claim.
	results := types.ValidationResults{
		OverallPass:     true,
		ComplianceLevel: "none",
		Tests:           []types.ValidationResult{},
		RegulatoryReady: false,
		AcademicReady:   false,
	}

	// Run configured statistical tests
	for _, testName := range f.config.StatisticalTests {
		result := f.runStatisticalTest(testName, data)
		results.Tests = append(results.Tests, result)
		if !result.Passed {
			results.OverallPass = false
		}
	}

	// Run population genetics tests
	for _, testName := range f.config.PopulationGeneticsTests {
		result := f.runPopulationGeneticsTest(testName, data)
		results.Tests = append(results.Tests, result)
		if !result.Passed {
			results.OverallPass = false
		}
	}

	// Run pipeline compatibility tests
	for _, testName := range f.config.PipelineCompatibilityTests {
		result := f.runPipelineCompatibilityTest(testName, data)
		results.Tests = append(results.Tests, result)
		if !result.Passed {
			results.OverallPass = false
		}
	}

	// Update summary and readiness based on results
	if len(results.Tests) == 0 {
		results.OverallPass = false
		results.Summary = "No validation tests configured"
	} else if results.OverallPass {
		results.Summary = fmt.Sprintf("All %d validation tests passed", len(results.Tests))
		results.AcademicReady = f.config.AcademicValidation.PeerReviewReady
		if results.AcademicReady {
			results.ComplianceLevel = "academic"
		} else {
			results.ComplianceLevel = "tests_passed"
		}
	} else {
		failedTests := 0
		for _, test := range results.Tests {
			if !test.Passed {
				failedTests++
			}
		}
		results.Summary = fmt.Sprintf("%d of %d tests failed", failedTests, len(results.Tests))
	}

	rc := f.config.RegulatoryCompliance
	if rc.FDACompliance || rc.EMACompliance || len(rc.ValidationStudies) > 0 || len(rc.QualityMetrics) > 0 {
		results.Summary += "; regulatory readiness not assessed (no FDA/EMA validation procedure is implemented)"
	}

	return results
}

// runStatisticalTest runs format_compliance, quality_distribution,
// sequence_complexity or read_length_consistency.
func (f *Framework) runStatisticalTest(testName string, data interface{}) types.ValidationResult {
	timestamp := time.Now()

	switch testName {
	case "format_compliance":
		return f.validateFASTQFormat(data, timestamp)
	case "quality_distribution":
		return f.validateQualityDistribution(data, timestamp)
	case "sequence_complexity":
		return f.validateSequenceComplexity(data, timestamp)
	case "read_length_consistency":
		return f.validateReadLengthConsistency(data, timestamp)
	default:
		return types.ValidationResult{
			TestName:     testName,
			Passed:       false,
			Confidence:   0.0,
			Message:      fmt.Sprintf("Unknown statistical test: %s", testName),
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}
}

// runPopulationGeneticsTest runs hardy_weinberg, base_composition or
// gc_content.
func (f *Framework) runPopulationGeneticsTest(testName string, data interface{}) types.ValidationResult {
	timestamp := time.Now()

	switch testName {
	case "hardy_weinberg":
		return f.validateHardyWeinberg(data, timestamp)
	case "base_composition":
		return f.validateBaseComposition(data, timestamp)
	case "gc_content":
		return f.validateGCContent(data, timestamp)
	default:
		return types.ValidationResult{
			TestName:     testName,
			Passed:       false,
			Confidence:   0.0,
			Message:      fmt.Sprintf("Unknown population genetics test: %s", testName),
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}
}

// runPipelineCompatibilityTest performs pipeline compatibility validation
func (f *Framework) runPipelineCompatibilityTest(testName string, data interface{}) types.ValidationResult {
	timestamp := time.Now()

	switch testName {
	case "fastq_parsing":
		return f.validateFASTQParsing(data, timestamp)
	case "quality_encoding":
		return f.validateQualityEncoding(data, timestamp)
	case "header_format":
		return f.validateHeaderFormat(data, timestamp)
	default:
		return types.ValidationResult{
			TestName:     testName,
			Passed:       false,
			Confidence:   0.0,
			Message:      fmt.Sprintf("Unknown pipeline compatibility test: %s", testName),
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}
}

// validateFASTQFormat checks only that the line count is a multiple of 4.
func (f *Framework) validateFASTQFormat(data interface{}, timestamp time.Time) types.ValidationResult {
	dataStr, ok := data.(string)
	if !ok {
		return types.ValidationResult{
			TestName:     "format_compliance",
			Passed:       false,
			Confidence:   0.0,
			Message:      "Invalid data type for FASTQ format validation",
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	lines := strings.Split(strings.TrimSpace(dataStr), "\n")
	if len(lines)%4 != 0 {
		return types.ValidationResult{
			TestName:     "format_compliance",
			Passed:       false,
			Confidence:   0.8,
			Message:      "FASTQ format violation: line count not divisible by 4",
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	return types.ValidationResult{
		TestName:     "format_compliance",
		Passed:       true,
		Confidence:   0.95,
		Message:      "FASTQ line count is a multiple of 4",
		Timestamp:    timestamp,
		ComplianceOK: true,
	}
}

// validateQualityDistribution checks that the mean Phred quality score is in a
// reasonable range (10–45) and that the distribution is not degenerate
// (all identical scores).
func (f *Framework) validateQualityDistribution(data interface{}, timestamp time.Time) types.ValidationResult {
	fastq := extractFASTQString(data)
	if fastq == "" {
		return types.ValidationResult{
			TestName:     "quality_distribution",
			Passed:       true,
			Confidence:   0.5,
			Message:      "Quality distribution validation skipped — no FASTQ data provided",
			Timestamp:    timestamp,
			ComplianceOK: true,
		}
	}

	lines := strings.Split(strings.TrimSpace(fastq), "\n")
	var scores []float64
	for i := 3; i < len(lines); i += 4 { // quality lines are at index 3, 7, 11, ...
		for _, ch := range lines[i] {
			scores = append(scores, float64(ch-33)) // Phred+33 decoding
		}
	}
	if len(scores) == 0 {
		return types.ValidationResult{
			TestName:     "quality_distribution",
			Passed:       true,
			Confidence:   0.5,
			Message:      "Quality distribution validation skipped — no quality lines found",
			Timestamp:    timestamp,
			ComplianceOK: true,
		}
	}

	sum := 0.0
	for _, s := range scores {
		sum += s
	}
	mean := sum / float64(len(scores))

	// Check for degenerate distribution (all scores identical)
	allSame := true
	for _, s := range scores {
		if s != scores[0] {
			allSame = false
			break
		}
	}
	if allSame {
		return types.ValidationResult{
			TestName:     "quality_distribution",
			Passed:       false,
			Confidence:   0.9,
			Message:      fmt.Sprintf("Degenerate quality distribution: all scores identical (%.0f)", scores[0]),
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	const minMeanPhred = 10.0
	const maxMeanPhred = 45.0
	if mean < minMeanPhred || mean > maxMeanPhred {
		return types.ValidationResult{
			TestName:     "quality_distribution",
			Passed:       false,
			Confidence:   0.9,
			Message:      fmt.Sprintf("Mean Phred quality %.1f outside expected range [%.0f, %.0f]", mean, minMeanPhred, maxMeanPhred),
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	return types.ValidationResult{
		TestName:     "quality_distribution",
		Passed:       true,
		Confidence:   0.85,
		Message:      fmt.Sprintf("Quality score distribution valid — mean Phred %.1f", mean),
		Timestamp:    timestamp,
		ComplianceOK: true,
	}
}

// validateSequenceComplexity fails if any read is a single repeated base.
func (f *Framework) validateSequenceComplexity(data interface{}, timestamp time.Time) types.ValidationResult {
	dataStr, ok := data.(string)
	if !ok {
		return types.ValidationResult{
			TestName:     "sequence_complexity",
			Passed:       false,
			Confidence:   0.0,
			Message:      "Invalid data type for sequence complexity validation",
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	// Basic complexity check: ensure sequences aren't all the same base
	lines := strings.Split(strings.TrimSpace(dataStr), "\n")
	for i := 1; i < len(lines); i += 4 { // Sequence lines
		if i < len(lines) {
			sequence := lines[i]
			if len(sequence) > 0 {
				firstBase := sequence[0]
				allSame := true
				for j := 0; j < len(sequence); j++ {
					if sequence[j] != firstBase {
						allSame = false
						break
					}
				}
				if allSame {
					return types.ValidationResult{
						TestName:     "sequence_complexity",
						Passed:       false,
						Confidence:   0.9,
						Message:      "Low sequence complexity detected: homopolymer sequence",
						Timestamp:    timestamp,
						ComplianceOK: false,
					}
				}
			}
		}
	}

	return types.ValidationResult{
		TestName:     "sequence_complexity",
		Passed:       true,
		Confidence:   0.8,
		Message:      "Sequence complexity validation passed",
		Timestamp:    timestamp,
		ComplianceOK: true,
	}
}

// validateReadLengthConsistency checks read lengths against f.readLengths.
func (f *Framework) validateReadLengthConsistency(data interface{}, timestamp time.Time) types.ValidationResult {
	result := func(passed bool, format string, args ...any) types.ValidationResult {
		return types.ValidationResult{
			TestName:     "read_length_consistency",
			Passed:       passed,
			Confidence:   0.9,
			Message:      fmt.Sprintf(format, args...),
			Timestamp:    timestamp,
			ComplianceOK: passed,
		}
	}
	fastq := extractFASTQString(data)
	if fastq == "" {
		return result(true, "Read length consistency validated")
	}

	lines := strings.Split(strings.TrimSpace(fastq), "\n")
	var lengths []int
	for i := 1; i < len(lines); i += 4 {
		lengths = append(lengths, len(strings.TrimSuffix(lines[i], "\r")))
	}
	if len(lengths) == 0 {
		return result(true, "Read length consistency validated")
	}
	st := summarizeLengths(lengths)
	exp := f.readLengths

	switch exp.Mode {
	case ReadLengthsAny:
		if st.min == st.max {
			return result(true, "All %d reads are %d bases", st.n, st.min)
		}
		return result(true, "Variable read lengths (informational; no expected distribution): %s", st)

	case ReadLengthsDistributed:
		name := exp.Name
		if name == "" {
			name = "expected"
		}
		var below, above int
		for _, l := range lengths {
			if exp.Min > 0 && l < exp.Min {
				below++
			}
			if exp.Max > 0 && l > exp.Max {
				above++
			}
		}
		if below+above > 0 {
			return result(false, "%d reads shorter than min %d and %d longer than max %d of the %s distribution: %s",
				below, exp.Min, above, exp.Max, name, st)
		}
		if exp.Mean > 0 && exp.SD > 0 {
			tol := 4 * exp.SD / math.Sqrt(float64(st.n))
			if math.Abs(st.mean-exp.Mean) > tol {
				return result(false, "Mean read length %.1f differs from the %s distribution mean %.1f by more than 4σ (%.1f): %s",
					st.mean, name, exp.Mean, tol, st)
			}
			return result(true, "Read lengths match the %s distribution (mean %.1f, expected %.1f ± %.1f): %s",
				name, st.mean, exp.Mean, tol, st)
		}
		return result(true, "Read lengths within the bounds of the %s distribution: %s", name, st)

	default: // ReadLengthsUniform
		expectedLen := lengths[0]
		inconsistent := 0
		for _, l := range lengths[1:] {
			if l != expectedLen {
				inconsistent++
			}
		}
		if inconsistent > 0 {
			return result(false, "Read length inconsistency: %d reads differ from expected length %d", inconsistent, expectedLen)
		}
		if exp.Length > 0 && expectedLen != exp.Length {
			return result(false, "All %d reads are %d bases, expected %d", st.n, expectedLen, exp.Length)
		}
		return result(true, "Read length consistency validated")
	}
}

// lengthStats summarizes read lengths.
type lengthStats struct {
	n, min, max, n50 int
	mean             float64
}

func (s lengthStats) String() string {
	return fmt.Sprintf("%d reads, min %d, max %d, mean %.1f, N50 %d", s.n, s.min, s.max, s.mean, s.n50)
}

// summarizeLengths returns count, bounds, mean and N50 (the length L such
// that reads of at least L bases hold half of all bases) of lengths.
func summarizeLengths(lengths []int) lengthStats {
	sorted := slices.Clone(lengths)
	slices.Sort(sorted)
	st := lengthStats{n: len(sorted), min: sorted[0], max: sorted[len(sorted)-1]}
	var total int64
	for _, l := range sorted {
		total += int64(l)
	}
	st.mean = float64(total) / float64(st.n)
	var acc int64
	for i := len(sorted) - 1; i >= 0; i-- {
		acc += int64(sorted[i])
		if 2*acc >= total {
			st.n50 = sorted[i]
			break
		}
	}
	return st
}

// validateHardyWeinberg runs a Hardy-Weinberg chi-squared test on the
// genotype counts in a *ValidationData. Without them it derives
// pseudo-genotypes from per-read GC content (computeGenotypeProxy); those are
// not genotypes, so that result does not test Hardy-Weinberg equilibrium at
// any real locus, and the message says so.
func (f *Framework) validateHardyWeinberg(data interface{}, timestamp time.Time) types.ValidationResult {
	var nHomRef, nHet, nHomAlt int
	var hasData bool
	source := "GC-content pseudo-genotypes, not real genotypes"

	switch v := data.(type) {
	case *ValidationData:
		if len(v.GenotypeCounts) > 0 {
			for _, gc := range v.GenotypeCounts {
				nHomRef += gc.NHomRef
				nHet += gc.NHet
				nHomAlt += gc.NHomAlt
			}
			hasData = true
			source = "genotype counts"
		} else if v.FASTQ != "" {
			nHomRef, nHet, nHomAlt, hasData = computeGenotypeProxy(v.FASTQ)
		}
	case string:
		if v != "" {
			nHomRef, nHet, nHomAlt, hasData = computeGenotypeProxy(v)
		}
	case []byte:
		if len(v) > 0 {
			nHomRef, nHet, nHomAlt, hasData = computeGenotypeProxy(string(v))
		}
	}

	if !hasData {
		return types.ValidationResult{
			TestName:     "hardy_weinberg",
			Passed:       true,
			Confidence:   0.5,
			Message:      "Hardy-Weinberg test skipped: no sequence data provided",
			Timestamp:    timestamp,
			ComplianceOK: true,
		}
	}

	result, err := statistics.TestHWE(nHomRef, nHet, nHomAlt)
	if err != nil {
		return types.ValidationResult{
			TestName:     "hardy_weinberg",
			Passed:       false,
			Confidence:   0.7,
			Message:      fmt.Sprintf("Hardy-Weinberg test error: %v", err),
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	message := fmt.Sprintf("Hardy-Weinberg χ²=%.3f p=%.4f (HomRef=%d, Het=%d, HomAlt=%d; %s)",
		result.ChiSquared, result.PValue, nHomRef, nHet, nHomAlt, source)

	return types.ValidationResult{
		TestName:     "hardy_weinberg",
		Passed:       result.PassesHWE,
		Confidence:   result.PValue,
		Message:      message,
		Timestamp:    timestamp,
		ComplianceOK: result.PassesHWE,
	}
}

// extractFASTQString extracts a FASTQ string from string, []byte, or *ValidationData.
// Returns "" when data is nil or an unsupported type.
func extractFASTQString(data interface{}) string {
	switch v := data.(type) {
	case string:
		return v
	case []byte:
		return string(v)
	case *ValidationData:
		return v.FASTQ
	}
	return ""
}

// computeGenotypeProxy derives pseudo-genotype counts from per-read GC fraction.
// Reads with GC more than 0.15 below the mean → HomRef, within ±0.15 → Het,
// above mean+0.15 → HomAlt. Requires at least 3 reads; returns ok=false otherwise.
func computeGenotypeProxy(fastq string) (nHomRef, nHet, nHomAlt int, ok bool) {
	lines := strings.Split(strings.TrimSpace(fastq), "\n")
	var gcFractions []float64

	for i := 1; i < len(lines); i += 4 {
		seq := strings.ToUpper(lines[i])
		if len(seq) == 0 {
			continue
		}
		gc, total := 0, 0
		for _, b := range seq {
			if b == 'G' || b == 'C' {
				gc++
			}
			if b == 'A' || b == 'T' || b == 'G' || b == 'C' {
				total++
			}
		}
		if total > 0 {
			gcFractions = append(gcFractions, float64(gc)/float64(total))
		}
	}

	if len(gcFractions) < 3 {
		return 0, 0, 0, false
	}

	mean := 0.0
	for _, f := range gcFractions {
		mean += f
	}
	mean /= float64(len(gcFractions))

	const threshold = 0.15
	for _, f := range gcFractions {
		switch {
		case f < mean-threshold:
			nHomRef++
		case f > mean+threshold:
			nHomAlt++
		default:
			nHet++
		}
	}

	return nHomRef, nHet, nHomAlt, true
}

// validateBaseComposition fails if any read contains a base other than
// A, C, G, T or N.
func (f *Framework) validateBaseComposition(data interface{}, timestamp time.Time) types.ValidationResult {
	dataStr, ok := data.(string)
	if !ok {
		return types.ValidationResult{
			TestName:     "base_composition",
			Passed:       false,
			Confidence:   0.0,
			Message:      "Invalid data type for base composition validation",
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	// Basic base composition check: ensure all bases are valid DNA
	lines := strings.Split(strings.TrimSpace(dataStr), "\n")
	for i := 1; i < len(lines); i += 4 { // Sequence lines
		if i < len(lines) {
			sequence := strings.ToUpper(lines[i])
			for _, base := range sequence {
				if base != 'A' && base != 'T' && base != 'G' && base != 'C' && base != 'N' {
					return types.ValidationResult{
						TestName:     "base_composition",
						Passed:       false,
						Confidence:   0.95,
						Message:      fmt.Sprintf("Invalid DNA base detected: %c", base),
						Timestamp:    timestamp,
						ComplianceOK: false,
					}
				}
			}
		}
	}

	return types.ValidationResult{
		TestName:     "base_composition",
		Passed:       true,
		Confidence:   0.95,
		Message:      "Base composition validation passed",
		Timestamp:    timestamp,
		ComplianceOK: true,
	}
}

// validateGCContent requires the overall GC fraction to be between 0.2 and
// 0.7.
func (f *Framework) validateGCContent(data interface{}, timestamp time.Time) types.ValidationResult {
	dataStr, ok := data.(string)
	if !ok {
		return types.ValidationResult{
			TestName:     "gc_content",
			Passed:       false,
			Confidence:   0.0,
			Message:      "Invalid data type for GC content validation",
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	// Calculate GC content
	lines := strings.Split(strings.TrimSpace(dataStr), "\n")
	totalBases := 0
	gcBases := 0

	for i := 1; i < len(lines); i += 4 { // Sequence lines
		if i < len(lines) {
			sequence := strings.ToUpper(lines[i])
			for _, base := range sequence {
				if base == 'G' || base == 'C' {
					gcBases++
				}
				if base == 'A' || base == 'T' || base == 'G' || base == 'C' {
					totalBases++
				}
			}
		}
	}

	if totalBases == 0 {
		return types.ValidationResult{
			TestName:     "gc_content",
			Passed:       false,
			Confidence:   0.5,
			Message:      "No valid DNA bases found for GC content calculation",
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	gcContent := float64(gcBases) / float64(totalBases)

	// Human genome GC content is about 41%; accept a wide 20-70% band.
	if gcContent < 0.2 || gcContent > 0.7 {
		return types.ValidationResult{
			TestName:     "gc_content",
			Passed:       false,
			Confidence:   0.8,
			Message:      fmt.Sprintf("GC content %.1f%% outside expected range (20-70%%)", gcContent*100),
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	return types.ValidationResult{
		TestName:     "gc_content",
		Passed:       true,
		Confidence:   0.9,
		Message:      fmt.Sprintf("GC content %.1f%% within expected range", gcContent*100),
		Timestamp:    timestamp,
		ComplianceOK: true,
	}
}

// validateFASTQParsing checks that the FASTQ data has valid 4-line record
// structure and that sequence and quality lines are the same length.
func (f *Framework) validateFASTQParsing(data interface{}, timestamp time.Time) types.ValidationResult {
	fastq := extractFASTQString(data)
	if fastq == "" {
		return types.ValidationResult{
			TestName:     "fastq_parsing",
			Passed:       true,
			Confidence:   0.5,
			Message:      "FASTQ parsing validation skipped — no FASTQ data provided",
			Timestamp:    timestamp,
			ComplianceOK: true,
		}
	}

	lines := strings.Split(strings.TrimSpace(fastq), "\n")
	if len(lines)%4 != 0 {
		return types.ValidationResult{
			TestName:     "fastq_parsing",
			Passed:       false,
			Confidence:   0.95,
			Message:      fmt.Sprintf("Invalid FASTQ structure: %d lines (not a multiple of 4)", len(lines)),
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	for i := 0; i+3 < len(lines); i += 4 {
		header, seq, sep, qual := lines[i], lines[i+1], lines[i+2], lines[i+3]
		if len(header) == 0 || header[0] != '@' {
			return types.ValidationResult{
				TestName:     "fastq_parsing",
				Passed:       false,
				Confidence:   0.95,
				Message:      fmt.Sprintf("Record %d: header line does not start with '@'", i/4+1),
				Timestamp:    timestamp,
				ComplianceOK: false,
			}
		}
		if len(sep) == 0 || sep[0] != '+' {
			return types.ValidationResult{
				TestName:     "fastq_parsing",
				Passed:       false,
				Confidence:   0.95,
				Message:      fmt.Sprintf("Record %d: separator line does not start with '+'", i/4+1),
				Timestamp:    timestamp,
				ComplianceOK: false,
			}
		}
		if len(seq) != len(qual) {
			return types.ValidationResult{
				TestName:     "fastq_parsing",
				Passed:       false,
				Confidence:   0.95,
				Message:      fmt.Sprintf("Record %d: sequence length %d ≠ quality length %d", i/4+1, len(seq), len(qual)),
				Timestamp:    timestamp,
				ComplianceOK: false,
			}
		}
	}

	return types.ValidationResult{
		TestName:     "fastq_parsing",
		Passed:       true,
		Confidence:   0.95,
		Message:      fmt.Sprintf("FASTQ parsing validated — %d records", len(lines)/4),
		Timestamp:    timestamp,
		ComplianceOK: true,
	}
}

// validateQualityEncoding verifies that all quality characters are in the
// Phred+33 ASCII range (33–126, i.e. '!' to '~').
func (f *Framework) validateQualityEncoding(data interface{}, timestamp time.Time) types.ValidationResult {
	fastq := extractFASTQString(data)
	if fastq == "" {
		return types.ValidationResult{
			TestName:     "quality_encoding",
			Passed:       true,
			Confidence:   0.5,
			Message:      "Quality encoding validation skipped — no FASTQ data provided",
			Timestamp:    timestamp,
			ComplianceOK: true,
		}
	}

	lines := strings.Split(strings.TrimSpace(fastq), "\n")
	for i := 3; i < len(lines); i += 4 {
		for pos, ch := range lines[i] {
			if ch < 33 || ch > 126 {
				return types.ValidationResult{
					TestName:     "quality_encoding",
					Passed:       false,
					Confidence:   0.95,
					Message:      fmt.Sprintf("Quality character out of Phred+33 range at record %d position %d: ASCII %d", i/4+1, pos, ch),
					Timestamp:    timestamp,
					ComplianceOK: false,
				}
			}
		}
	}

	return types.ValidationResult{
		TestName:     "quality_encoding",
		Passed:       true,
		Confidence:   0.95,
		Message:      "Quality encoding (Phred+33) validated — all characters in ASCII 33–126",
		Timestamp:    timestamp,
		ComplianceOK: true,
	}
}

// validateHeaderFormat checks only that every header line starts with '@'.
func (f *Framework) validateHeaderFormat(data interface{}, timestamp time.Time) types.ValidationResult {
	dataStr, ok := data.(string)
	if !ok {
		return types.ValidationResult{
			TestName:     "header_format",
			Passed:       false,
			Confidence:   0.0,
			Message:      "Invalid data type for header format validation",
			Timestamp:    timestamp,
			ComplianceOK: false,
		}
	}

	lines := strings.Split(strings.TrimSpace(dataStr), "\n")
	for i := 0; i < len(lines); i += 4 { // Header lines
		if i < len(lines) {
			header := lines[i]
			if !strings.HasPrefix(header, "@") {
				return types.ValidationResult{
					TestName:     "header_format",
					Passed:       false,
					Confidence:   0.95,
					Message:      "Invalid FASTQ header format: missing '@' prefix",
					Timestamp:    timestamp,
					ComplianceOK: false,
				}
			}
		}
	}

	return types.ValidationResult{
		TestName:     "header_format",
		Passed:       true,
		Confidence:   0.95,
		Message:      "Header format validation passed",
		Timestamp:    timestamp,
		ComplianceOK: true,
	}
}
