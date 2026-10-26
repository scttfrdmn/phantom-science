package cmd

import (
	"fmt"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/pkg/types"
	"github.com/spf13/cobra"
)

// Fingerprinting settings of the scenario, signal and transcriptome
// commands (#31) come from the same places as generate's: the sequencer
// config (--config: fingerprinting, enterprise.customer_id, dataset_id), a
// phantom config (-p: enterprise.fingerprinting and customer_id, which take
// precedence) and --dataset-id.

// applyFingerprintFlags resolves cfg's fingerprinting settings from the
// phantom config (-p) and --dataset-id.
func applyFingerprintFlags(cmd *cobra.Command, cfg *types.SequencerConfig) error {
	var phantom *types.PhantomConfig
	if path, _ := cmd.InheritedFlags().GetString("phantom"); path != "" {
		p, err := cli.NewConfigLoader(false).LoadPhantomConfig(path)
		if err != nil {
			return fmt.Errorf("loading phantom config: %w", err)
		}
		phantom = p
	}
	cfg.ResolveFingerprinting(phantom)
	if cmd.Flags().Changed("dataset-id") {
		cfg.DatasetID, _ = cmd.Flags().GetString("dataset-id")
	}
	return nil
}

// commandRunSeed resolves the run seed of a command without a sequencer
// configuration of its own (signal, transcriptome) from --seed and the
// fingerprinting settings, and logs a drawn seed or salt to logw.
func commandRunSeed(cmd *cobra.Command, seed uint64, logf func(format string, args ...any)) (genomics.RunSeed, error) {
	cfg := &types.SequencerConfig{}
	if path, _ := cmd.InheritedFlags().GetString("config"); path != "" {
		loaded, err := cli.NewConfigLoader(false).LoadSequencerConfig(path)
		if err != nil {
			return genomics.RunSeed{}, fmt.Errorf("loading config: %w", err)
		}
		cfg = loaded
	}
	if err := applyFingerprintFlags(cmd, cfg); err != nil {
		return genomics.RunSeed{}, err
	}
	cfg.Seed = seed
	rs, err := genomics.ResolveRunSeed(*cfg, cmd.Name())
	if err != nil {
		return genomics.RunSeed{}, err
	}
	switch {
	case rs.Drawn && rs.Fingerprint != nil:
		logf("   [info] fingerprint salt %d (keyed seed; pass --seed %d with the same key, customer and dataset to reproduce this run)\n", rs.Fingerprint.Salt, rs.Fingerprint.Salt)
	case rs.Drawn:
		logf("   [info] seed %d (pass --seed %d to reproduce this run)\n", rs.Seed, rs.Seed)
	}
	return rs, nil
}
