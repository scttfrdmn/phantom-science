package fingerprinting

import "strings"

// Read-ID fingerprint tags. Generation writes them as SAM-style optional
// fields in the FASTQ header comment, after the platform comment and any
// RG:Z: read group, separated by tabs:
//
//	@<read name> 1:N:0:1	RG:Z:<sample>	XC:Z:<customer>	XD:Z:<dataset>
//
// so the read name and the Illumina comment stay intact, and aligners that
// copy FASTQ comments to SAM (bwa mem -C, minimap2 -y) carry the tags into
// the alignment records. X? tags are reserved for local use by the SAM
// specification.
const (
	CustomerTag = "XC" // customer identifier
	DatasetTag  = "XD" // dataset identifier
)

// ReadTags returns the read-ID fingerprint fields for customerID and
// datasetID, e.g. ["XC:Z:acme", "XD:Z:ds1"]; an empty identifier gives no
// field. Whitespace in a value (which would split the FASTQ comment) is
// replaced by '_'.
func ReadTags(customerID, datasetID string) []string {
	var tags []string
	if customerID != "" {
		tags = append(tags, CustomerTag+":Z:"+tagValue(customerID))
	}
	if datasetID != "" {
		tags = append(tags, DatasetTag+":Z:"+tagValue(datasetID))
	}
	return tags
}

// ParseReadTags returns the customer and dataset of the read-ID fingerprint
// fields in a FASTQ header line (with or without the leading '@'); missing
// tags give empty strings.
func ParseReadTags(header string) (customerID, datasetID string) {
	for _, f := range strings.Split(header, "\t")[1:] {
		switch {
		case strings.HasPrefix(f, CustomerTag+":Z:"):
			customerID = f[len(CustomerTag)+3:]
		case strings.HasPrefix(f, DatasetTag+":Z:"):
			datasetID = f[len(DatasetTag)+3:]
		}
	}
	return customerID, datasetID
}

// tagValue makes v a single SAM Z value without whitespace.
func tagValue(v string) string {
	return strings.Map(func(r rune) rune {
		switch r {
		case ' ', '\t', '\n', '\r', '\v', '\f':
			return '_'
		}
		return r
	}, v)
}
