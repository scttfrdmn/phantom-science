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
	key, err := fingerprinting.ResolveKey(c.Fingerprinting.KeyFile)
	switch {
	case err != nil:
		e.fpErr = fmt.Errorf("fingerprinting: %w", err)
		return
	case key == nil:
		e.fpErr = fmt.Errorf("fingerprinting needs a key: set fingerprinting.key_file or %s (a keyed seed without a key proves nothing); "+
			"to write only the read-ID labels, set steganography_layers: [%q]", fingerprinting.KeyEnvVar, fingerprinting.LayerReadIDs)
		return
	case c.Enterprise.CustomerID == "":
		e.fpErr = fmt.Errorf("fingerprinting needs a customer ID: set enterprise.customer_id")
		return
	case c.DatasetID == "":
		e.fpErr = fmt.Errorf("fingerprinting needs a dataset ID: set dataset_id (or generate --dataset-id)")
		return
	}
	seed, err := fingerprinting.DeriveSeed(key, c.Enterprise.CustomerID, c.DatasetID, salt)
	if err != nil {
		e.fpErr = fmt.Errorf("fingerprinting: %w", err)
		return
	}
	e.seed = seed
	e.fp = keyedSeed{enabled: true, customer: c.Enterprise.CustomerID, dataset: c.DatasetID, salt: salt, saltDrawn: drawn}
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
