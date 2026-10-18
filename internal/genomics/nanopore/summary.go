package nanopore

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"path/filepath"
	"strconv"
	"sync"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

// SummaryFileName is MinKNOW's default sequencing summary file name.
const SummaryFileName = "sequencing_summary.txt"

// SummaryColumns are the sequencing_summary.txt columns written, a subset of
// MinKNOW's that covers what PycoQC needs.
var SummaryColumns = []string{
	"filename_fastq", "read_id", "run_id", "channel", "start_time", "duration",
	"sequence_length_template", "mean_qscore_template",
}

// SummaryWriter writes one sequencing_summary.txt row per observed read. It
// implements the engine's ReadObserver and is safe for concurrent use.
//
// start_time is in seconds since the run start and duration in seconds,
// both with microsecond precision, matching the read header's start_time.
// mean_qscore_template is MeanQScore of the qualities as written to FASTQ
// (capped at maxPhred).
type SummaryWriter struct {
	info     RunInfo
	maxPhred byte
	open     func(rec *seqread.Record) (io.WriteCloser, error)
	fileName func(file string) string

	mu     sync.Mutex
	opened bool
	wc     io.WriteCloser
	bw     *bufio.Writer
	err    error
}

// NewSummaryWriter returns a writer that opens its output lazily, on the
// first read, by calling open with that read (so the destination can depend
// on whether reads go to files). open may return a nil writer to write no
// summary. fileName maps a read's output file to the file name written in
// filename_fastq (e.g. adding ".gz" for files compressed afterwards); nil
// uses the base name.
func NewSummaryWriter(info RunInfo, maxPhred byte, open func(rec *seqread.Record) (io.WriteCloser, error), fileName func(file string) string) *SummaryWriter {
	if fileName == nil {
		fileName = func(f string) string {
			if f == "" {
				return ""
			}
			return filepath.Base(f)
		}
	}
	return &SummaryWriter{info: info, maxPhred: maxPhred, open: open, fileName: fileName}
}

// ObserveRead appends rec's row.
func (s *SummaryWriter) ObserveRead(rec *seqread.Record) error {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.err != nil {
		return s.err
	}
	if !s.opened {
		s.opened = true
		wc, err := s.open(rec)
		if err != nil {
			s.err = fmt.Errorf("opening sequencing summary: %w", err)
			return s.err
		}
		if wc == nil {
			return nil
		}
		s.wc, s.bw = wc, bufio.NewWriterSize(wc, 256*1024)
		for i, c := range SummaryColumns {
			if i > 0 {
				s.bw.WriteByte('\t') // #nosec G104 -- bufio errors are sticky and surface on Flush
			}
			s.bw.WriteString(c) // #nosec G104
		}
		s.bw.WriteByte('\n') // #nosec G104
	}
	if s.bw == nil {
		return nil
	}
	var buf [64]byte
	b := s.bw
	b.WriteString(s.fileName(rec.File)) // #nosec G104 -- sticky; checked on Flush
	b.WriteByte('\t')
	b.WriteString(rec.Name)
	b.WriteByte('\t')
	b.WriteString(s.info.RunID)
	b.WriteByte('\t')
	b.Write(strconv.AppendInt(buf[:0], int64(rec.Channel), 10))
	b.WriteByte('\t')
	b.Write(appendSeconds(buf[:0], rec.StartTime.Sub(s.info.Start)))
	b.WriteByte('\t')
	b.Write(appendSeconds(buf[:0], rec.Duration))
	b.WriteByte('\t')
	b.Write(strconv.AppendInt(buf[:0], int64(len(rec.Seq)), 10))
	b.WriteByte('\t')
	b.Write(strconv.AppendFloat(buf[:0], MeanQScore(rec.Qual, s.maxPhred), 'f', 6, 64))
	if err := b.WriteByte('\n'); err != nil {
		s.err = fmt.Errorf("writing sequencing summary: %w", err)
		return s.err
	}
	return nil
}

// Close flushes and closes the summary.
func (s *SummaryWriter) Close() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.wc == nil {
		return s.err
	}
	err := s.bw.Flush()
	err = errors.Join(err, s.wc.Close())
	s.wc, s.bw = nil, nil
	if err != nil {
		s.err = errors.Join(s.err, fmt.Errorf("writing sequencing summary: %w", err))
	}
	return s.err
}

// appendSeconds formats d as seconds with exactly six decimals, without
// floating-point rounding.
func appendSeconds(b []byte, d time.Duration) []byte {
	us := int64(d / time.Microsecond)
	if us < 0 {
		b = append(b, '-')
		us = -us
	}
	b = strconv.AppendInt(b, us/1e6, 10)
	b = append(b, '.')
	frac := us % 1e6
	for div := int64(1e5); div > 0; div /= 10 {
		b = append(b, byte('0'+frac/div%10))
	}
	return b
}
