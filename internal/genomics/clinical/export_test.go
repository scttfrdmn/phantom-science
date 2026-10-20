//go:build !release

// export_test.go exposes internal functions for black-box tests in package clinical_test.
package clinical

import (
	"io"
	mathrand "math/rand"

	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// ExportGenerateHaplotypeFASTQ exposes generateHaplotypeFASTQ for testing.
func ExportGenerateHaplotypeFASTQ(
	w io.Writer,
	engine *genomics.Engine,
	cfg *types.SequencerConfig,
	trio *TrioHaplotypes,
	member string,
	region *genomics.GenomicRegion,
	readLength int,
	readCount int,
) (int, error) {
	return generateHaplotypeFASTQ(w, engine, cfg, trio, member, region, readLength, readCount)
}

// ExportSimulateCrossover exposes simulateCrossover for testing.
func ExportSimulateCrossover(h1, h2 []int, positions []int, geneticMapCM float64, rng *mathrand.Rand) []int {
	return simulateCrossover(h1, h2, positions, geneticMapCM, rng)
}
