package genomics

// Regeneration for keyed-seed verification (#27).
//
// A unit's reads are a pure function of (flow-cell seed, configuration,
// batch, mate, number of reads), and a unit cut to k reads is the first k
// reads of the full unit (see the file plan in pipeline.go). Regenerate
// therefore reproduces any file of a run without writing it, one file at a
// time, and hands each read to a callback; memory does not grow with the
// run.

import (
	"errors"
	"fmt"
	"io"

	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

// RunFile identifies one FASTQ file of a run.
type RunFile struct {
	FlowCell int // 1-based
	Batch    int // 1-based unit (batch) number
	Mate     int // 1, or 2 for the second file of a pair
	Reads    int // reads (read pairs) in the file
}

// RegenOptions selects what Regenerate reproduces.
type RegenOptions struct {
	// Files lists the files to regenerate, in order (for example from a
	// run manifest). Nil regenerates the configured plan: every flow cell's
	// units, as Generate writes them. With a total_bases target, units
	// are regenerated until their bases reach the target, which may be one
	// unit more than Generate wrote (a superset is harmless for matching).
	Files []RunFile
	// MaxFiles stops after this many files; 0 for no limit.
	MaxFiles int
	// MaxReads stops after this many reads (each mate counted); 0 for no
	// limit.
	MaxReads int64
}

// errRegenDone stops regeneration at a limit.
var errRegenDone = errors.New("regeneration limit reached")

// regenObserver passes every written read to a callback.
type regenObserver struct {
	file    RunFile
	observe func(RunFile, []byte) error
	reads   *int64
	bases   *int64
	max     int64
}

func (o *regenObserver) ObserveRead(rec *seqread.Record) error {
	*o.reads++
	*o.bases += int64(len(rec.Seq))
	if err := o.observe(o.file, rec.Seq); err != nil {
		return err
	}
	if o.max > 0 && *o.reads >= o.max {
		return errRegenDone
	}
	return nil
}

func (o *regenObserver) Close() error { return nil }

// Regenerate reproduces the run's FASTQ files without writing them and calls
// observe with the bases of every read, in file order. It returns the number
// of files regenerated (a file stopped by MaxReads counts). The engine must
// be configured like the run (seed included); Regenerate writes nothing and
// needs no output directory.
func (e *Engine) Regenerate(opts RegenOptions, observe func(f RunFile, seq []byte) error) (int, error) {
	if err := e.checkExtensions(); err != nil {
		return 0, err
	}
	n := max(e.config.FlowCells, 1)
	engines := make([]*Engine, 0, n)
	seen := map[string]bool{}
	for k := 1; k <= n; k++ {
		fc, _, err := e.newFlowCellEngine(k, n, seen)
		if err != nil {
			return 0, err
		}
		engines = append(engines, fc)
	}

	var reads int64
	files := 0
	// one regenerates f; done reports a limit reached.
	one := func(f RunFile) (done bool, bases int64, err error) {
		if f.FlowCell < 1 || f.FlowCell > n {
			return false, 0, fmt.Errorf("flow cell %d: the run has %d", f.FlowCell, n)
		}
		if f.Mate != 1 && !(f.Mate == 2 && e.config.PairedEnd) {
			return false, 0, fmt.Errorf("mate %d: not a file of this run", f.Mate)
		}
		if f.Batch < 1 || f.Reads < 1 {
			return false, 0, fmt.Errorf("batch %d with %d reads: not a file of this run", f.Batch, f.Reads)
		}
		fc := engines[f.FlowCell-1]
		obs := &regenObserver{file: f, observe: observe, reads: &reads, bases: &bases, max: opts.MaxReads}
		job := fileJob{path: fc.outputPath(f.Batch, f.Mate), batch: f.Batch, mate: f.Mate, reads: f.Reads, obs: obs, plain: true}
		err = fc.writeFileStream(io.Discard, nil, job)
		files++
		if errors.Is(err, errRegenDone) {
			return true, bases, nil
		}
		if err != nil {
			return false, bases, fmt.Errorf("regenerating flow cell %d batch %d mate %d: %w", f.FlowCell, f.Batch, f.Mate, err)
		}
		return opts.MaxFiles > 0 && files >= opts.MaxFiles, bases, nil
	}

	if opts.Files != nil {
		for _, f := range opts.Files {
			if done, _, err := one(f); err != nil || done {
				return files, err
			}
		}
		return files, nil
	}
	for k, fc := range engines {
		p := fc.outputPlan()
		var total int64
		for unit := 1; p.has(unit); unit++ {
			for _, mate := range fc.mates() {
				done, bases, err := one(RunFile{FlowCell: k + 1, Batch: unit, Mate: mate, Reads: p.reads(unit)})
				if err != nil || done {
					return files, err
				}
				total += bases
			}
			if (p.bases > 0 && total >= p.bases) || (p.cap > 0 && total >= p.cap) {
				break
			}
		}
	}
	return files, nil
}
