package genomics

// Keyed-seed fingerprints (docs/CHARTER.md, principle 2; #27).
//
// With fingerprinting enabled and the keyed_seed layer (the default when no
// layer is listed), the run seed is not SequencerConfig.Seed but
//
//	fingerprinting.DeriveSeed(key, enterprise.customer_id, dataset_id, salt)
//
// where salt is SequencerConfig.Seed (--seed) or, when that is zero, a
// random value drawn at run start and logged. Nothing else changes: every
// value of the output follows from the derived seed through the usual
// seeded streams, so the output is byte-identical to an unfingerprinted run
// with seed = derived seed, for any number of workers. The manifest records
// the scheme, customer, dataset and salt, never the derived seed or the key.
//
// Fingerprinting fails before any output is written without a key, without
// a customer or dataset ID, with an unknown layer, or with the
// quality_scores layer, which edited quality values and is no longer
// written. The read_ids layer alone writes only the read-ID labels, with no
// keyed seed.

import (
	"fmt"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// keyedSeed describes the run's keyed-seed fingerprint.
type keyedSeed struct {
	enabled           bool
	customer, dataset string
	salt              uint64
	saltDrawn         bool // the salt was drawn, not configured; it is logged
}

// resolveSeed sets the run seed from the configuration: the configured or a
// drawn seed, or, with a keyed-seed fingerprint, the seed derived from the
// key. Unusable fingerprint settings are recorded in fpErr, which stops
// generation before any output is written.
func (e *Engine) resolveSeed() {
	c := e.config
	e.seed, e.seedDrawn, e.seedErr = ResolveSeed(c.Seed)
	if !c.Fingerprinting.Enable {
		return
	}
	keyed, err := fingerprinting.Layers(c.Fingerprinting.SteganographyLayers)
	if err != nil {
		e.fpErr = err
		return
	}
	if !keyed {
		return
	}
	if e.seedErr != nil {
		return
	}
	salt, drawn := e.seed, e.seedDrawn
	e.seedDrawn = false // the derived seed is never logged
	seed, err := deriveKeyedSeed(c, salt)
	if err != nil {
		e.fpErr = err
		return
	}
	e.seed = seed
	e.fp = keyedSeed{enabled: true, customer: c.Enterprise.CustomerID, dataset: c.DatasetID, salt: salt, saltDrawn: drawn}
}

// deriveKeyedSeed returns the keyed-seed run seed of configuration c for
// salt: the key from c.Fingerprinting.KeyFile or the environment, the
// customer from c.Enterprise.CustomerID and the dataset from c.DatasetID,
// all required.
func deriveKeyedSeed(c types.SequencerConfig, salt uint64) (uint64, error) {
	key, err := fingerprinting.ResolveKey(c.Fingerprinting.KeyFile)
	switch {
	case err != nil:
		return 0, fmt.Errorf("fingerprinting: %w", err)
	case key == nil:
		return 0, fmt.Errorf("fingerprinting needs a key: set fingerprinting.key_file or %s (a keyed seed without a key proves nothing); "+
			"to write only the read-ID labels, set steganography_layers: [%q]", fingerprinting.KeyEnvVar, fingerprinting.LayerReadIDs)
	case c.Enterprise.CustomerID == "":
		return 0, fmt.Errorf("fingerprinting needs a customer ID: set enterprise.customer_id")
	case c.DatasetID == "":
		return 0, fmt.Errorf("fingerprinting needs a dataset ID: set dataset_id (or generate --dataset-id)")
	}
	seed, err := fingerprinting.DeriveSeed(key, c.Enterprise.CustomerID, c.DatasetID, salt)
	if err != nil {
		return 0, fmt.Errorf("fingerprinting: %w", err)
	}
	return seed, nil
}

// RunSeed is the seed of a run of a command other than generate (scenario,
// signal, transcriptome), resolved by ResolveRunSeed.
type RunSeed struct {
	// Seed is the run seed: the configured or drawn seed, or, with a
	// keyed-seed fingerprint, the seed derived from the key. A derived seed
	// must not be logged or recorded.
	Seed uint64
	// Drawn reports that the seed (or, when fingerprinted, the salt) was
	// drawn rather than configured; it should be logged so the run can be
	// reproduced: Seed when Fingerprint is nil, else Fingerprint.Salt.
	Drawn bool
	// Fingerprint is the manifest's fingerprint record; nil without a
	// keyed-seed fingerprint.
	Fingerprint *manifest.Fingerprint
}

// ResolveRunSeed resolves the run seed of command (for messages) from c:
// c.Seed (zero draws one) and, when c.Fingerprinting is enabled, the
// keyed-seed fingerprint, with c.Seed (or the drawn value) as the salt.
// Every random value of such a run must come from the returned seed, so a
// fingerprinted run equals an unfingerprinted run with Seed = the derived
// seed. The commands write no read headers with labels, so the read_ids
// layer and the customer and dataset label settings are errors.
func ResolveRunSeed(c types.SequencerConfig, command string) (RunSeed, error) {
	seed, drawn, err := ResolveSeed(c.Seed)
	if err != nil {
		return RunSeed{}, err
	}
	if !c.Fingerprinting.Enable {
		return RunSeed{Seed: seed, Drawn: drawn}, nil
	}
	keyed, err := fingerprinting.Layers(c.Fingerprinting.SteganographyLayers)
	if err != nil {
		return RunSeed{}, err
	}
	if !keyed || c.Fingerprinting.CustomerFingerprint || c.Fingerprinting.DatasetFingerprint {
		return RunSeed{}, fmt.Errorf("%s supports only the %q fingerprint layer, without read-ID labels (customer_fingerprint, dataset_fingerprint, the %q layer)",
			command, fingerprinting.LayerKeyedSeed, fingerprinting.LayerReadIDs)
	}
	derived, err := deriveKeyedSeed(c, seed)
	if err != nil {
		return RunSeed{}, err
	}
	return RunSeed{Seed: derived, Drawn: drawn, Fingerprint: &manifest.Fingerprint{
		Scheme:   fingerprinting.KeyedSeedScheme,
		Customer: c.Enterprise.CustomerID,
		Dataset:  c.DatasetID,
		Salt:     seed,
	}}, nil
}

// Fingerprinted reports whether the run's seed is a keyed-seed fingerprint,
// and its salt.
func (e *Engine) Fingerprinted() (salt uint64, ok bool) {
	return e.fp.salt, e.fp.enabled
}

// manifestFingerprint returns the manifest's fingerprint record, or nil.
func (e *Engine) manifestFingerprint() *manifest.Fingerprint {
	if !e.fp.enabled {
		return nil
	}
	return &manifest.Fingerprint{
		Scheme:   fingerprinting.KeyedSeedScheme,
		Customer: e.fp.customer,
		Dataset:  e.fp.dataset,
		Salt:     e.fp.salt,
	}
}

// recordedSeed returns seed as the manifest and object tags may show it:
// a keyed-seed run's seeds (the run's and those derived from it) are not
// recorded, since anyone could otherwise regenerate the run without the key.
func (e *Engine) recordedSeed(seed uint64) uint64 {
	if e.fp.enabled || (e.parent != nil && e.parent.fp.enabled) {
		return 0
	}
	return seed
}
