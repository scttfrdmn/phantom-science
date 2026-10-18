package fingerprinting

import (
	"fmt"
	"strconv"
	"strings"
	"time"
)

// Fingerprint is the parsed form of a fingerprint payload,
// "<customer>:<dataset>:stegano-<version>:ts-<unix seconds>".
type Fingerprint struct {
	Customer       string
	Dataset        string
	SteganoVersion string
	Created        time.Time
}

// FormatFingerprint returns the payload string for a fingerprint.
func FormatFingerprint(customer, dataset, steganoVersion string, created time.Time) string {
	return fmt.Sprintf("%s:%s:stegano-%s:ts-%d", customer, dataset, steganoVersion, created.Unix())
}

// ParseFingerprint parses a payload written by FormatFingerprint. The
// customer is the first field and the dataset is everything up to the
// version field, so a dataset ID may contain ':' but a customer ID may not.
func ParseFingerprint(s string) (Fingerprint, error) {
	parts := strings.Split(s, ":")
	if len(parts) < 4 {
		return Fingerprint{}, fmt.Errorf("fingerprint %q: want customer:dataset:stegano-<version>:ts-<unix>", s)
	}
	n := len(parts)
	version, ok := strings.CutPrefix(parts[n-2], "stegano-")
	if !ok {
		return Fingerprint{}, fmt.Errorf("fingerprint %q: missing stegano-<version> field", s)
	}
	tsField, ok := strings.CutPrefix(parts[n-1], "ts-")
	if !ok {
		return Fingerprint{}, fmt.Errorf("fingerprint %q: missing ts-<unix> field", s)
	}
	ts, err := strconv.ParseInt(tsField, 10, 64)
	if err != nil {
		return Fingerprint{}, fmt.Errorf("fingerprint %q: bad timestamp: %w", s, err)
	}
	return Fingerprint{
		Customer:       parts[0],
		Dataset:        strings.Join(parts[1:n-2], ":"),
		SteganoVersion: version,
		Created:        time.Unix(ts, 0).UTC(),
	}, nil
}
