//go:build !release

package cli

import (
	"context"

	"github.com/scttfrdmn/phantom-science/internal/genomics/clinical"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// scenarioRegenerate regenerates a scenario run (clinical.Regenerate). The
// scenario command, and so its regeneration, is not in release builds.
func scenarioRegenerate(ctx context.Context, cfg types.SequencerConfig, sc *manifest.Scenario, observe func(member string, seq []byte) error) (int, error) {
	return clinical.Regenerate(ctx, cfg, sc, observe)
}
