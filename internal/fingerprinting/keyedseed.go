package fingerprinting

import (
	"crypto/hmac"
	"crypto/sha256"
	"encoding/binary"
	"errors"
	"fmt"
	"strings"
)

// Keyed-seed fingerprints (docs/CHARTER.md, principle 2).
//
// A fingerprinted run does not edit its output. Instead its seed is derived
// from the fingerprint key, the customer, the dataset and a salt:
//
//	seed = first 8 bytes, big-endian, of
//	       HMAC-SHA256(key, "phantom-science/fingerprint/v1\x00" ‖ customer ‖
//	                        "\x00" ‖ dataset ‖ "\x00" ‖ salt)
//
// where salt is encoded as 8 bytes, big-endian. Every base, quality, length
// and header value follows from that seed exactly as from any other seed, so
// a fingerprinted run is byte-identical to an unfingerprinted run with the
// derived seed, and its statistics are those of any run. Provenance is
// proven by regenerating the run from the key, the configuration, the
// customer, the dataset and the salt, and matching suspect reads against it
// (see Matcher). Without the key the seed cannot be computed, so the data
// cannot be attributed, and nothing in the data can be stripped.

// Fingerprinting layers (FingerprintConfig.SteganographyLayers).
const (
	// LayerKeyedSeed derives the run seed from the key (the default when
	// fingerprinting is enabled without layers).
	LayerKeyedSeed = "keyed_seed"
	// LayerReadIDs names the read-ID labels (customer_fingerprint,
	// dataset_fingerprint) as the only layer: the labels are written
	// without a keyed seed. The labels are visible header text, not a
	// fingerprint.
	LayerReadIDs = "read_ids"
	// LayerQualityScores is the quality-score embedding of phantom-science
	// 0.9.0 and earlier. It edits quality values, so generation rejects it;
	// verify still extracts it from old files.
	LayerQualityScores = "quality_scores"
)

// KeyedSeedScheme identifies the keyed-seed derivation in run manifests.
const KeyedSeedScheme = "keyed_seed/v1"

// keyedSeedDomain separates this derivation from any other use of the key.
const keyedSeedDomain = "phantom-science/fingerprint/v1\x00"

// ErrQualityScoresLayer rejects the quality_scores layer at generation.
var ErrQualityScoresLayer = errors.New(`steganography layer "quality_scores" is no longer supported for generation: it changes quality values, ` +
	`which a fingerprint must never do (docs/CHARTER.md, principle 2). Use "keyed_seed" (the default when fingerprinting is enabled), ` +
	`which leaves the data unchanged; "phantom-science verify" still reads quality-score fingerprints from files made by phantom-science 0.9.0 and earlier`)

// Layers reports what the configured layers enable: keyed is true for the
// keyed-seed fingerprint, which is the default for an empty list. The
// quality_scores layer and unknown layer names are errors.
func Layers(layers []string) (keyed bool, err error) {
	if len(layers) == 0 {
		return true, nil
	}
	for _, l := range layers {
		switch strings.TrimSpace(l) {
		case LayerKeyedSeed:
			keyed = true
		case LayerReadIDs:
		case LayerQualityScores:
			return false, ErrQualityScoresLayer
		default:
			return false, fmt.Errorf("unknown steganography layer %q (known layers: %q, %q)", l, LayerKeyedSeed, LayerReadIDs)
		}
	}
	return keyed, nil
}

// DeriveSeed returns the keyed-seed run seed for key, customer, dataset and
// salt. customer and dataset must be non-empty and contain no NUL byte, and
// key must pass ValidateKey. A derived value of zero (probability 2^-64) is
// returned as 1, because a zero seed means "draw one".
func DeriveSeed(key []byte, customer, dataset string, salt uint64) (uint64, error) {
	if err := ValidateKey(key); err != nil {
		return 0, err
	}
	for name, v := range map[string]string{"customer ID": customer, "dataset ID": dataset} {
		if v == "" {
			return 0, fmt.Errorf("keyed-seed fingerprint: %s is empty", name)
		}
		if strings.IndexByte(v, 0) >= 0 {
			return 0, fmt.Errorf("keyed-seed fingerprint: %s contains a NUL byte", name)
		}
	}
	mac := hmac.New(sha256.New, key)
	mac.Write([]byte(keyedSeedDomain)) // #nosec G104 -- hash writes cannot fail
	mac.Write([]byte(customer))        // #nosec G104
	mac.Write([]byte{0})               // #nosec G104
	mac.Write([]byte(dataset))         // #nosec G104
	mac.Write([]byte{0})               // #nosec G104
	var s [8]byte
	binary.BigEndian.PutUint64(s[:], salt)
	mac.Write(s[:]) // #nosec G104
	seed := binary.BigEndian.Uint64(mac.Sum(nil)[:8])
	if seed == 0 {
		seed = 1
	}
	return seed, nil
}
