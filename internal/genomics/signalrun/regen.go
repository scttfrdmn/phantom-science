// Package signalrun holds the parts of the signal command (raw nanopore
// signal, BLOW5) that verification shares with it: reading the input reads,
// mapping read IDs, loading the profile and pore model, and regenerating
// the signal of a run (#31).
//
// The raw signal of a read is a pure function of the run seed, the read ID
// and the sequence it is synthesised from (signal.SynthesizeRead), so a
// signal run is regenerated read by read from its inputs, profile, pore
// model and seed, in any order and with any number of workers. A keyed-seed
// fingerprint is the run seed: a fingerprinted signal run is the
// unfingerprinted run whose seed is the derived seed.
package signalrun

import (
	"errors"
	"fmt"
	"io"
	"runtime"
	"sync"

	"github.com/scttfrdmn/phantom-science/internal/genomics/signal"
)

// LoadModel returns the profile named profile with the overrides sets
// applied, and the pore model at modelPath read in units ("auto", "pa" or
// "standardised"; auto detects them), with the units used. A k that
// differs from the profile's is reported to logw.
func LoadModel(profile string, sets map[string]float64, modelPath, units string, logw io.Writer) (signal.Profile, *signal.PoreModel, string, error) {
	p, err := signal.Preset(profile)
	if err != nil {
		return signal.Profile{}, nil, "", err
	}
	for name, v := range sets {
		if err := p.Set(name, v); err != nil {
			return signal.Profile{}, nil, "", err
		}
	}
	if err := p.Validate(); err != nil {
		return signal.Profile{}, nil, "", err
	}
	if units == "auto" {
		if units, err = DetectModelUnits(modelPath); err != nil {
			return signal.Profile{}, nil, "", err
		}
	}
	var lopts []signal.LoadOption
	if units == "standardised" {
		lopts = append(lopts, signal.WithScaling(p.ModelScale, p.ModelShift), signal.WithDefaultStdv(p.ModelStdv))
	}
	model, err := signal.LoadPoreModel(modelPath, lopts...)
	if err != nil {
		return signal.Profile{}, nil, "", err
	}
	if p.ModelK > 0 && model.K != p.ModelK {
		fmt.Fprintf(logw, "   [warn] pore model k=%d, profile %s expects k=%d\n", model.K, p.Name, p.ModelK)
	}
	return p, model, units, nil
}

// RegenSpec describes a signal run to regenerate: the paths of its input
// reads (in order), its truth file ("" for none), pore model and the
// parameters its manifest records.
type RegenSpec struct {
	Inputs     []string
	TruthPath  string
	ModelPath  string
	ModelUnits string // "pa" or "standardised" (as recorded)
	Profile    string
	Sets       map[string]float64
	// Workers synthesise signal in parallel (0: number of CPUs); the
	// signals do not depend on it, only the order observe sees them in.
	Workers int
	// MaxReads stops after this many reads (0: all).
	MaxReads int64
}

// ErrStop, returned by an observe function, ends Regenerate early without
// an error.
var ErrStop = errors.New("stop regeneration")

// Regenerate synthesises the raw signal of every read of the run spec
// describes with seed, without writing it, and calls observe with each
// read's ID and raw samples (from one goroutine, in no particular order).
// It returns the number of reads synthesised.
func Regenerate(spec RegenSpec, seed uint64, observe func(id string, raw []int16) error) (int64, error) {
	p, model, _, err := LoadModel(spec.Profile, spec.Sets, spec.ModelPath, spec.ModelUnits, io.Discard)
	if err != nil {
		return 0, err
	}
	var truth map[string][]byte
	if spec.TruthPath != "" {
		if truth, _, err = LoadTruth(spec.TruthPath); err != nil {
			return 0, err
		}
	}
	workers := spec.Workers
	if workers <= 0 {
		workers = runtime.NumCPU()
	}
	type job struct {
		id  string
		seq []byte
	}
	type done struct {
		id  string
		raw []int16
		err error
	}
	jobs := make(chan job, 2*workers)
	results := make(chan done, 2*workers)
	quit := make(chan struct{})
	var wg sync.WaitGroup
	for range workers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobs {
				rd, err := signal.SynthesizeRead(seed, j.id, j.seq, p, model)
				var d done
				if err != nil {
					d = done{id: j.id, err: fmt.Errorf("read %s: %w", j.id, err)}
				} else {
					d = done{id: j.id, raw: rd.Raw}
				}
				select {
				case results <- d:
				case <-quit:
				}
			}
		}()
	}
	var readErr error
	go func() {
		defer close(jobs)
		var queued int64
		for _, path := range spec.Inputs {
			s, err := OpenSeqFile(path)
			if err != nil {
				readErr = err
				return
			}
			for {
				name, _, seq, err := s.Next()
				if err == io.EOF {
					break
				}
				if err != nil {
					s.Close()
					readErr = err
					return
				}
				id, _, src, _ := Source(name, seq, truth)
				if len(src) < model.K {
					continue
				}
				select {
				case jobs <- job{id, src}:
				case <-quit:
					s.Close()
					return
				}
				if queued++; spec.MaxReads > 0 && queued >= spec.MaxReads {
					s.Close()
					return
				}
			}
			s.Close()
		}
	}()
	go func() {
		wg.Wait()
		close(results)
	}()
	var n int64
	var err2 error
	stopped := false
	for d := range results {
		if stopped {
			continue
		}
		switch {
		case d.err != nil:
			err2 = d.err
		default:
			n++
			err2 = observe(d.id, d.raw)
		}
		if err2 != nil {
			stopped = true
			close(quit)
		}
	}
	if errors.Is(err2, ErrStop) {
		err2 = nil
	}
	if stopped {
		return n, err2
	}
	return n, readErr
}
