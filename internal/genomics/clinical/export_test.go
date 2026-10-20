//go:build !release

// export_test.go exposes internal functions for black-box tests in package clinical_test.
package clinical

import (
	mathrand "math/rand"

	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// ExportMemberHaplotypes returns member's two haplotype sequences built from
// the window bases base (0-based start) and trio, as a trio scenario does
// without a reference (the simulated REF alleles are written into a copy of
// base first).
func ExportMemberHaplotypes(base []byte, start int64, trio *TrioHaplotypes, member string) ([2][]byte, error) {
	tp := exportPlan(base, start, trio, 0)
	return tp.memberHaplotypes(member)
}

// ExportWriteMemberReads writes member's reads, sampled by the engine with
// cfg (ReadsPerBatch reads or pairs) from its haplotypes as
// ExportMemberHaplotypes builds them, to r1 (and r2 for separate mates). It
// returns the number of FASTQ records.
func ExportWriteMemberReads(cfg types.SequencerConfig, base []byte, start int64, trio *TrioHaplotypes, member, r1, r2 string) (int, error) {
	cfg.TotalBatches = 1
	tp := exportPlan(base, start, trio, cfg.Seed)
	return writeMemberReads(cfg, tp, member, output{r1, r2})
}

func exportPlan(base []byte, start int64, trio *TrioHaplotypes, seed uint64) *trioPlan {
	b := append([]byte(nil), base...)
	writeSyntheticREF(trio, b, start)
	return &trioPlan{trio: trio, base: b, start: start, runSeed: seed}
}

// ExportSimulateCrossover exposes simulateCrossover for testing.
func ExportSimulateCrossover(h1, h2 []int, positions []int, geneticMapCM float64, rng *mathrand.Rand) []int {
	return simulateCrossover(h1, h2, positions, geneticMapCM, rng)
}
