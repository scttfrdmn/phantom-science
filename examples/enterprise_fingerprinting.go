// Command enterprise_fingerprinting generates a small fingerprinted FASTQ
// dataset, then proves where its files came from. It does in Go what these
// commands do:
//
//	phantom-science generate --config ... -p configs/phantom/v1.0/enterprise-settings.yaml
//	phantom-science verify --key-file fingerprint.key --manifest OUT/phantom-manifest.json FILE...
//
// The fingerprint is a keyed seed (docs/CHARTER.md, principle 2): the run
// seed is HMAC-SHA256(key, customer, dataset, salt), so the output is an
// ordinary run, byte for byte, that only the key holder can attribute. The
// data is not edited: no base, quality or header changes. Verification
// regenerates the run from the key and the manifest (customer, dataset,
// salt and configuration) and matches the files' reads against it.
//
// Usage:
//
//	openssl rand -base64 32 > fingerprint.key && chmod 600 fingerprint.key
//	go run ./examples -key-file fingerprint.key [-out DIR] [-customer ID] [-dataset ID] [-reads N]
//
// Without -key-file the key comes from $PHANTOM_FINGERPRINT_KEY. No stegano
// library or CGO is needed.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/scttfrdmn/phantom-science/internal/cli"
	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/internal/genomics/manifest"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Options configures FingerprintDataset.
type Options struct {
	OutputDir  string
	CustomerID string
	DatasetID  string
	KeyFile    string // empty: $PHANTOM_FINGERPRINT_KEY
	Reads      int    // read pairs in the single output file pair
	Seed       uint64 // the salt; 0: draw one (logged to stderr, recorded in the manifest)
}

// FingerprintDataset generates one paired-end FASTQ file pair in
// o.OutputDir, fingerprinted for o.CustomerID and o.DatasetID with a keyed
// seed, and returns the paths of the files and of the run manifest.
func FingerprintDataset(o Options) (files []string, manifestPath string, err error) {
	if err := os.MkdirAll(o.OutputDir, 0o750); err != nil {
		return nil, "", fmt.Errorf("creating output directory: %w", err)
	}
	config := types.SequencerConfig{
		SampleName:    "fingerprint-example",
		Lane:          1,
		OutputDir:     o.OutputDir,
		ReadLength:    150,
		ReadsPerBatch: o.Reads,
		TotalBatches:  1,
		PairedEnd:     true,
		BaseQuality:   35,
		QualityDecay:  0.1,
		DatasetID:     o.DatasetID,
		Seed:          o.Seed,
		Enterprise:    types.EnterpriseConfig{CustomerID: o.CustomerID},
		Fingerprinting: types.FingerprintConfig{
			Enable:              true,
			SteganographyLayers: []string{fingerprinting.LayerKeyedSeed},
			KeyFile:             o.KeyFile,
		},
	}
	engine := genomics.NewEngine(config)
	defer engine.Close()
	if err := engine.Generate(); err != nil {
		return nil, "", err
	}
	return engine.GeneratedFiles(), filepath.Join(o.OutputDir, manifest.FileName), nil
}

// VerifyDataset checks, with the key from keyFile (or
// $PHANTOM_FINGERPRINT_KEY), that every file's reads come from the run that
// manifestPath describes, reports one verdict per file to w, and returns an
// error unless every file matches.
func VerifyDataset(w io.Writer, files []string, keyFile, manifestPath string) error {
	return cli.VerifyFiles(w, files, cli.VerifyOptions{KeyFile: keyFile, Manifests: []string{manifestPath}})
}

func main() {
	out := flag.String("out", "./fingerprinted_output", "output directory")
	customer := flag.String("customer", "ACME_GENOMICS_001", "customer ID the run is fingerprinted for")
	dataset := flag.String("dataset", "cancer_study_001", "dataset ID the run is fingerprinted for")
	keyFile := flag.String("key-file", "", "fingerprint key file, mode 600 (default: $"+fingerprinting.KeyEnvVar+")")
	reads := flag.Int("reads", 1000, "read pairs to generate")
	flag.Parse()

	files, man, err := FingerprintDataset(Options{
		OutputDir: *out, CustomerID: *customer, DatasetID: *dataset,
		KeyFile: *keyFile, Reads: *reads,
	})
	if err != nil {
		fmt.Fprintf(os.Stderr, "generation failed: %v\n", err)
		os.Exit(1)
	}
	if err := VerifyDataset(os.Stdout, files, *keyFile, man); err != nil {
		fmt.Fprintf(os.Stderr, "verification failed: %v\n", err)
		os.Exit(1)
	}
}
