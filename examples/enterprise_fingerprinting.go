// Command enterprise_fingerprinting generates a small FASTQ dataset whose
// files carry an encrypted steganographic fingerprint, then extracts and
// checks the fingerprint of every file. It does in Go what these commands do:
//
//	phantom-science generate --config ... -p configs/phantom/v1.0/enterprise-settings.yaml
//	phantom-science verify --customer ID --dataset ID FILE...
//
// Each output file is fingerprinted once, after it is written, with stegano's
// file-level embedding into the quality scores. The fingerprint payload is
// "<customer>:<dataset>:stegano-<version>:ts-<unix time>", and it is
// encrypted with a key (stegano's Argon2id key derivation and AES-256-GCM),
// so it can only be read with the same key.
//
// Usage:
//
//	openssl rand -base64 32 > fingerprint.key && chmod 600 fingerprint.key
//	go run ./examples -key-file fingerprint.key [-out DIR] [-customer ID] [-dataset ID] [-reads N]
//
// Without -key-file the key comes from $PHANTOM_FINGERPRINT_KEY. Embedding
// and extraction need the stegano library: build with CGO after
// `make stegano` (and CGO_LDFLAGS as set by the Makefile). In builds without
// it, generation fails before any file is written.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/fingerprinting"
	"github.com/scttfrdmn/phantom-science/internal/genomics"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// Options configures FingerprintDataset.
type Options struct {
	OutputDir  string
	CustomerID string
	DatasetID  string
	KeyFile    string // empty: $PHANTOM_FINGERPRINT_KEY
	Reads      int    // read pairs in the single output file pair
	Seed       uint64 // 0: draw a seed (logged to stderr)
}

// FingerprintDataset generates one paired-end FASTQ file pair in
// o.OutputDir with an encrypted quality-score fingerprint naming
// o.CustomerID and o.DatasetID, and returns the paths of the files written.
func FingerprintDataset(o Options) ([]string, error) {
	if err := os.MkdirAll(o.OutputDir, 0o750); err != nil {
		return nil, fmt.Errorf("creating output directory: %w", err)
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
			SteganographyLayers: []string{"quality_scores"},
			KeyFile:             o.KeyFile,
		},
	}
	engine := genomics.NewEngine(config)
	defer engine.Close()
	if err := engine.Generate(); err != nil {
		return nil, err
	}
	return engine.GeneratedFiles(), nil
}

// VerifyDataset extracts the fingerprint of each file with the key from
// keyFile (or $PHANTOM_FINGERPRINT_KEY), reports it to w, and returns an
// error unless every file has an encrypted fingerprint naming customerID and
// datasetID.
func VerifyDataset(w io.Writer, files []string, keyFile, customerID, datasetID string) error {
	key, err := fingerprinting.ResolveKey(keyFile)
	if err != nil {
		return err
	}
	var opts []fingerprinting.Option
	if key != nil {
		opts = append(opts, fingerprinting.WithKey(key))
	}
	verifier, err := fingerprinting.NewVerifier(opts...)
	if err != nil {
		return err
	}
	failed := 0
	for _, path := range files {
		problem := ""
		res, err := verifier.Verify(path)
		switch {
		case err != nil:
			problem = err.Error()
		case !res.Found:
			problem = "no fingerprint"
		case !res.Encrypted:
			problem = "fingerprint is not encrypted"
		default:
			fp, perr := fingerprinting.ParseFingerprint(res.Fingerprint)
			switch {
			case perr != nil:
				problem = perr.Error()
			case fp.Customer != customerID || fp.Dataset != datasetID:
				problem = fmt.Sprintf("fingerprint names customer %q dataset %q", fp.Customer, fp.Dataset)
			default:
				fmt.Fprintf(w, "ok   %s: customer=%s dataset=%s created=%s stegano=%s (encrypted)\n",
					filepath.Base(path), fp.Customer, fp.Dataset, fp.Created.UTC().Format(time.RFC3339), fp.SteganoVersion)
				continue
			}
		}
		fmt.Fprintf(w, "FAIL %s: %s\n", filepath.Base(path), problem)
		failed++
	}
	if failed > 0 {
		return fmt.Errorf("%d of %d files failed verification", failed, len(files))
	}
	return nil
}

func main() {
	out := flag.String("out", "./fingerprinted_output", "output directory")
	customer := flag.String("customer", "ACME_GENOMICS_001", "customer ID written into the fingerprint")
	dataset := flag.String("dataset", "cancer_study_001", "dataset ID written into the fingerprint")
	keyFile := flag.String("key-file", "", "encryption key file, mode 600 (default: $"+fingerprinting.KeyEnvVar+")")
	reads := flag.Int("reads", 1000, "read pairs to generate")
	flag.Parse()

	files, err := FingerprintDataset(Options{
		OutputDir: *out, CustomerID: *customer, DatasetID: *dataset,
		KeyFile: *keyFile, Reads: *reads,
	})
	if err != nil {
		fmt.Fprintf(os.Stderr, "generation failed: %v\n", err)
		os.Exit(1)
	}
	if err := VerifyDataset(os.Stdout, files, *keyFile, *customer, *dataset); err != nil {
		fmt.Fprintf(os.Stderr, "verification failed: %v\n", err)
		os.Exit(1)
	}
}
