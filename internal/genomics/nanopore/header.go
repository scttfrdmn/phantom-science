package nanopore

import (
	"math"
	"math/rand/v2"
	"strconv"
	"strings"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

// Channel timing. Durations are whole microseconds, so header start times
// (printed to the microsecond) and summary offsets agree exactly.
const (
	minGap     = 100 * time.Millisecond  // pore idle time before a read (capture)
	maxGap     = 1500 * time.Millisecond //
	minAdapter = 100 * time.Millisecond  // adapter / leader time added to each read
	maxAdapter = 300 * time.Millisecond  //
)

// StartTimeLayout formats header start times like MinKNOW: RFC 3339 with
// microseconds and a numeric UTC offset ("2024-03-01T12:00:00.250000+00:00").
const StartTimeLayout = "2006-01-02T15:04:05.000000-07:00"

// Placement is where and when one read was sequenced.
type Placement struct {
	Channel  int           // 1..Channels
	Read     int           // per-channel read number (MinKNOW read=)
	Start    time.Duration // offset of the read start from the run start
	Duration time.Duration // read duration
}

// ChannelModel assigns the reads of one output stream to channels and gives
// each a per-channel read number, start time and duration.
//
// Reads are dealt to channels in rounds: each round of Channels reads visits
// every channel once, in a fresh random order. On a channel, a read starts
// after an idle gap following the previous read and lasts
// length/SpeedBPS plus a short adapter time, so reads on a channel never
// overlap and start in increasing order.
//
// Streams are independent (each draws only from its own generator), so they
// can be generated in any order. Stream s (0-based) owns the time window
// [s·W, (s+1)·W) and the per-channel read numbers (s·K, (s+1)·K], where
// K = ceil(readsPerStream/Channels) is the number of rounds and W = K times
// the longest slot (maximum gap + maximum adapter time + the duration of
// slotBases bases). Every read of up to slotBases bases therefore ends
// inside its stream's window, and per-channel read numbers and start times
// increase across the whole run. A read longer than slotBases stays
// non-overlapping within its stream but may run into the next stream's
// window on the same channel.
type ChannelModel struct {
	info      RunInfo
	offset    time.Duration // start of this stream's window
	firstRead int           // read numbers of this stream start after this
	clocks    []time.Duration
	perm      []int
	n         int
}

// NewChannelModel returns the channel model for stream (0-based) of a run
// whose streams hold readsPerStream reads, sizing the time windows for reads
// of up to slotBases bases. info must have Channels and SpeedBPS set.
func NewChannelModel(info RunInfo, stream, readsPerStream, slotBases int) *ChannelModel {
	c := info.Channels
	rounds := (max(readsPerStream, 1) + c - 1) / c
	slot := maxGap + maxAdapter + basesDuration(slotBases, info.SpeedBPS)
	m := &ChannelModel{
		info:      info,
		offset:    time.Duration(stream) * time.Duration(rounds) * slot,
		firstRead: stream * rounds,
		clocks:    make([]time.Duration, c),
	}
	for i := range m.clocks {
		m.clocks[i] = m.offset
	}
	return m
}

// basesDuration is the translocation time of n bases, in whole microseconds.
func basesDuration(n int, speed float64) time.Duration {
	return time.Duration(math.Round(float64(n)/speed*1e6)) * time.Microsecond
}

// uniformMicros draws a duration uniformly from [lo, hi] in whole microseconds.
func uniformMicros(r *rand.Rand, lo, hi time.Duration) time.Duration {
	span := int64((hi - lo) / time.Microsecond)
	return lo + time.Duration(r.Int64N(span+1))*time.Microsecond
}

// Next places the next read, of length bases, drawing only from r.
func (m *ChannelModel) Next(r *rand.Rand, length int) Placement {
	c := len(m.clocks)
	round, idx := m.n/c, m.n%c
	if idx == 0 {
		m.perm = r.Perm(c)
	}
	m.n++
	ch := m.perm[idx]
	start := m.clocks[ch] + uniformMicros(r, minGap, maxGap)
	dur := basesDuration(length, m.info.SpeedBPS) + uniformMicros(r, minAdapter, maxAdapter)
	m.clocks[ch] = start + dur
	return Placement{Channel: ch + 1, Read: m.firstRead + round + 1, Start: start, Duration: dur}
}

// FormatComment returns the MinKNOW header fields after the read ID, in
// MinKNOW order:
//
//	runid=<hex> read=<n> ch=<channel> start_time=<RFC 3339> flow_cell_id=<id> protocol_group_id=<id> sample_id=<id>
func FormatComment(info RunInfo, p Placement) string {
	var b strings.Builder
	b.Grow(160 + len(info.RunID) + len(info.FlowCellID) + len(info.ProtocolGroupID) + len(info.SampleID))
	b.WriteString("runid=")
	b.WriteString(info.RunID)
	b.WriteString(" read=")
	b.WriteString(strconv.Itoa(p.Read))
	b.WriteString(" ch=")
	b.WriteString(strconv.Itoa(p.Channel))
	b.WriteString(" start_time=")
	b.WriteString(info.Start.Add(p.Start).Format(StartTimeLayout))
	b.WriteString(" flow_cell_id=")
	b.WriteString(info.FlowCellID)
	b.WriteString(" protocol_group_id=")
	b.WriteString(info.ProtocolGroupID)
	b.WriteString(" sample_id=")
	b.WriteString(info.SampleID)
	return b.String()
}

// phredError[q] is the error probability of Phred score q.
var phredError = func() (t [256]float64) {
	for q := range t {
		t[q] = math.Pow(10, -float64(q)/10)
	}
	return t
}()

// MeanQScore returns the read's mean quality as ONT basecallers define it:
// the Phred score of the mean per-base error probability,
// -10·log10(mean(10^(-q/10))). Scores above maxPhred count as maxPhred (the
// value written to FASTQ). An empty read scores 0.
func MeanQScore(qual []byte, maxPhred byte) float64 {
	if len(qual) == 0 {
		return 0
	}
	sum := 0.0
	for _, q := range qual {
		sum += phredError[min(q, maxPhred)]
	}
	return float64(-10 * math.Log10(sum/float64(len(qual))))
}

// NewHeaderFormatter returns a per-stream read-header function for the
// genomics pipeline: it draws a UUID4 read ID and the read's placement from
// r, records Channel, StartTime and Duration on rec (for the sequencing
// summary), and returns the read ID and MinKNOW header fields. The stream
// index is taken from the first read's 1-based Batch. Call it once per
// output stream.
func NewHeaderFormatter(info RunInfo, readsPerStream, slotBases int) func(r *rand.Rand, rec *seqread.Record) (name, comment string) {
	var m *ChannelModel
	return func(r *rand.Rand, rec *seqread.Record) (string, string) {
		if m == nil {
			m = NewChannelModel(info, max(rec.Batch-1, 0), readsPerStream, slotBases)
		}
		id := NewUUID4(r)
		p := m.Next(r, len(rec.Seq))
		rec.Channel, rec.StartTime, rec.Duration = p.Channel, info.Start.Add(p.Start), p.Duration
		return id, FormatComment(info, p)
	}
}
