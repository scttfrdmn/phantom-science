//go:build release

package cli

import (
	"context"
	"errors"

	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// scenarioRegenerate fails in release builds, which have no scenario
// command.
func scenarioRegenerate(context.Context, types.SequencerConfig, *manifest.Scenario, func(string, []byte) error) (int, error) {
	return 0, errors.New("a scenario run is verified with a build that has the scenario command (built without the release tag), as the run was made")
}
