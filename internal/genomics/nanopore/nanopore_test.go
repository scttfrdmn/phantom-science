//go:build !release

package nanopore

import (
	"bytes"
	"io"
	"math"
	"math/rand/v2"
	"regexp"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/seqread"
)

func testRand(seed uint64) *rand.Rand { return rand.New(rand.NewPCG(seed, 7)) }

// uuid4RE is a canonical lowercase UUID with version 4 and the RFC 4122
// variant (8, 9, a or b).
var uuid4RE = regexp.MustCompile(`^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$`)

// headerRE is a MinKNOW FASTQ header, with fields in MinKNOW order.
var headerRE = regexp.MustCompile(`^@(\S+) runid=([0-9a-f]{40}) read=(\d+) ch=(\d+) start_time=(\S+) flow_cell_id=(\S+) protocol_group_id=(\S+) sample_id=(\S+)$`)

func TestNewUUID4_VersionVariantAndUniqueness(t *testing.T) {
	r := testRand(1)
	seen := make(map[string]struct{}, 100000)
	for i := 0; i < 100000; i++ {
		u := NewUUID4(r)
		if !uuid4RE.MatchString(u) {
			t.Fatalf("UUID %q is not a canonical version-4 RFC 4122 UUID", u)
		}
		if _, dup := seen[u]; dup {
			t.Fatalf("duplicate UUID %q after %d draws", u, i)
		}
		seen[u] = struct{}{}
	}
	if NewUUID4(testRand(9)) != NewUUID4(testRand(9)) {
		t.Error("UUIDs are not reproducible from the generator")
	}
}

func TestRunInfo_DefaultsFromSeed(t *testing.T) {
	a, b := RunInfo{}.WithDefaults(42), RunInfo{}.WithDefaults(42)
	if a != b {
		t.Errorf("defaults differ for the same seed: %+v vs %+v", a, b)
	}
	if c := (RunInfo{}).WithDefaults(43); c.RunID == a.RunID || c.FlowCellID == a.FlowCellID && c.Start.Equal(a.Start) {
		t.Errorf("defaults do not depend on the seed: %+v vs %+v", a, c)
	}
	if !regexp.MustCompile(`^[0-9a-f]{40}$`).MatchString(a.RunID) {
		t.Errorf("run ID %q is not 40 hex characters", a.RunID)
	}
	if !regexp.MustCompile(`^PA[A-Z]\d{5}$`).MatchString(a.FlowCellID) {
		t.Errorf("flow cell ID %q is not a PromethION ID", a.FlowCellID)
	}
	if a.Channels != DefaultChannels || a.SpeedBPS != DefaultSpeedBPS || a.SampleID != DefaultSampleID || a.ProtocolGroupID == "" {
		t.Errorf("unexpected defaults %+v", a)
	}
	if err := a.Validate(); err != nil {
		t.Errorf("defaults do not validate: %v", err)
	}
	// Set fields are kept, and the generated ones do not change.
	set := RunInfo{RunID: "abc", Channels: 512}.WithDefaults(42)
	if set.RunID != "abc" || set.Channels != 512 || set.FlowCellID != a.FlowCellID || !set.Start.Equal(a.Start) {
		t.Errorf("WithDefaults overrode set fields or changed other defaults: %+v", set)
	}
	for _, bad := range []RunInfo{
		{RunID: "a b"}, {FlowCellID: "x=y"}, {SampleID: "s\t1"}, {Channels: -1}, {SpeedBPS: -400},
	} {
		if err := bad.WithDefaults(1).Validate(); err == nil {
			t.Errorf("Validate accepted %+v", bad)
		}
	}
}

// fakeRead runs one read through a header formatter like the pipeline does.
func fakeRead(f func(*rand.Rand, *seqread.Record) (string, string), r *rand.Rand, batch, n int) seqread.Record {
	rec := seqread.Record{Batch: batch, Seq: bytes.Repeat([]byte("A"), n), Qual: bytes.Repeat([]byte{20}, n)}
	rec.Name, rec.Comment = f(r, &rec)
	return rec
}

func TestHeaderFormatter_RoundTrip(t *testing.T) {
	info := RunInfo{SampleID: "S1", ProtocolGroupID: "exp_7"}.WithDefaults(5)
	f := NewHeaderFormatter(info, 100, 2000)
	r := testRand(2)
	for i := 0; i < 100; i++ {
		rec := fakeRead(f, r, 1, 1000+i)
		m := headerRE.FindStringSubmatch(rec.Header())
		if m == nil {
			t.Fatalf("header %q does not match the MinKNOW layout", rec.Header())
		}
		if !uuid4RE.MatchString(m[1]) || m[1] != rec.Name {
			t.Errorf("read ID %q", m[1])
		}
		if m[2] != info.RunID || m[6] != info.FlowCellID || m[7] != "exp_7" || m[8] != "S1" {
			t.Errorf("run fields %q", rec.Header())
		}
		if ch, _ := strconv.Atoi(m[4]); ch != rec.Channel || ch < 1 || ch > info.Channels {
			t.Errorf("ch=%s, record channel %d", m[4], rec.Channel)
		}
		st, err := time.Parse(time.RFC3339Nano, m[5])
		if err != nil {
			t.Fatalf("start_time %q is not RFC 3339: %v", m[5], err)
		}
		if !st.Equal(rec.StartTime) {
			t.Errorf("header start_time %v != record %v", st, rec.StartTime)
		}
		if want := time.Duration(math.Round(float64(len(rec.Seq))/info.SpeedBPS*1e6)) * time.Microsecond; rec.Duration < want+minAdapter || rec.Duration > want+maxAdapter {
			t.Errorf("duration %v for %d bases at %g b/s", rec.Duration, len(rec.Seq), info.SpeedBPS)
		}
	}
}

type span struct {
	read       int
	start, end time.Duration
}

// checkChannels fails if reads on a channel overlap, or if start times or
// read numbers do not increase on a channel.
func checkChannels(t *testing.T, byCh map[int][]span) {
	t.Helper()
	for ch, s := range byCh {
		for i := 1; i < len(s); i++ {
			if s[i].start < s[i-1].end {
				t.Fatalf("channel %d: read at %v starts before the previous read ends at %v", ch, s[i].start, s[i-1].end)
			}
			if s[i].read <= s[i-1].read {
				t.Fatalf("channel %d: read number %d after %d", ch, s[i].read, s[i-1].read)
			}
		}
	}
}

func TestChannelModel_NoOverlapWithinAndAcrossStreams(t *testing.T) {
	info := RunInfo{Channels: 7}.WithDefaults(3)
	const perStream, length = 30, 5000
	byCh := map[int][]span{}
	// Streams are generated out of order: windows must not depend on order.
	for _, stream := range []int{2, 0, 1} {
		m := NewChannelModel(info, stream, perStream, length)
		r := testRand(uint64(stream))
		for i := 0; i < perStream; i++ {
			p := m.Next(r, length)
			byCh[p.Channel] = append(byCh[p.Channel], span{p.Read, p.Start, p.Start + p.Duration})
		}
	}
	for ch, s := range byCh {
		// Sort by start; the check then verifies non-overlap across streams.
		for i := 1; i < len(s); i++ {
			for j := i; j > 0 && s[j].start < s[j-1].start; j-- {
				s[j], s[j-1] = s[j-1], s[j]
			}
		}
		byCh[ch] = s
		if len(s) < 3*(perStream/7) || len(s) > 3*((perStream+6)/7) {
			t.Errorf("channel %d got %d reads: rounds should spread reads evenly", ch, len(s))
		}
	}
	checkChannels(t, byCh)
}

func TestChannelModel_LongReadsStayOrderedWithinStream(t *testing.T) {
	info := RunInfo{Channels: 3}.WithDefaults(3)
	m := NewChannelModel(info, 0, 50, 100) // window sized for 100 bp reads
	r := testRand(4)
	byCh := map[int][]span{}
	for i := 0; i < 50; i++ {
		p := m.Next(r, 1000+r.IntN(100000))
		byCh[p.Channel] = append(byCh[p.Channel], span{p.Read, p.Start, p.Start + p.Duration})
	}
	checkChannels(t, byCh)
}

type nopCloser struct{ io.Writer }

func (nopCloser) Close() error { return nil }

func TestSummaryWriter_RowsMatchHeaders(t *testing.T) {
	info := RunInfo{SampleID: "s"}.WithDefaults(11)
	f := NewHeaderFormatter(info, 20, 300)
	var buf bytes.Buffer
	opens := 0
	sw := NewSummaryWriter(info, 50, func(*seqread.Record) (io.WriteCloser, error) {
		opens++
		return nopCloser{&buf}, nil
	}, nil)
	r := testRand(6)
	var recs []seqread.Record
	for i := 0; i < 20; i++ {
		rec := fakeRead(f, r, 1, 50+i*10)
		rec.File = "/out/x.fastq"
		for j := range rec.Qual {
			rec.Qual[j] = byte(5 + (i+j)%60) // includes scores above the cap
		}
		if err := sw.ObserveRead(&rec); err != nil {
			t.Fatal(err)
		}
		recs = append(recs, rec)
	}
	if err := sw.Close(); err != nil {
		t.Fatal(err)
	}
	if opens != 1 {
		t.Errorf("summary opened %d times", opens)
	}
	lines := strings.Split(strings.TrimSuffix(buf.String(), "\n"), "\n")
	if lines[0] != strings.Join(SummaryColumns, "\t") {
		t.Fatalf("header row %q", lines[0])
	}
	if len(lines) != 1+len(recs) {
		t.Fatalf("%d rows for %d reads", len(lines)-1, len(recs))
	}
	for i, rec := range recs {
		c := strings.Split(lines[i+1], "\t")
		if c[0] != "x.fastq" || c[1] != rec.Name || c[2] != info.RunID || c[3] != strconv.Itoa(rec.Channel) || c[6] != strconv.Itoa(len(rec.Seq)) {
			t.Errorf("row %d %q does not match read %s ch %d len %d", i, lines[i+1], rec.Name, rec.Channel, len(rec.Seq))
		}
		st, _ := strconv.ParseFloat(c[4], 64)
		if got := info.Start.Add(time.Duration(math.Round(st*1e6)) * time.Microsecond); !got.Equal(rec.StartTime) {
			t.Errorf("row %d start_time %s is %v, header says %v", i, c[4], got, rec.StartTime)
		}
		if d, _ := strconv.ParseFloat(c[5], 64); time.Duration(math.Round(d*1e6))*time.Microsecond != rec.Duration {
			t.Errorf("row %d duration %s, record %v", i, c[5], rec.Duration)
		}
		// Recompute the ONT mean qscore by hand from the capped qualities.
		var sum float64
		for _, q := range rec.Qual {
			sum += math.Pow(10, -float64(min(q, 50))/10)
		}
		want := -10 * math.Log10(sum/float64(len(rec.Qual)))
		if got, _ := strconv.ParseFloat(c[7], 64); math.Abs(got-want) > 1e-5 {
			t.Errorf("row %d mean_qscore %s, want %.6f", i, c[7], want)
		}
	}
}

func TestSummaryWriter_NilDestinationWritesNothing(t *testing.T) {
	sw := NewSummaryWriter(RunInfo{}.WithDefaults(1), 50, func(*seqread.Record) (io.WriteCloser, error) { return nil, nil }, nil)
	rec := seqread.Record{Name: "r", Seq: []byte("AC"), Qual: []byte{1, 2}}
	if err := sw.ObserveRead(&rec); err != nil {
		t.Fatal(err)
	}
	if err := sw.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestMeanQScore(t *testing.T) {
	// Q10 and Q30: error probabilities 0.1 and 0.001 average to 0.0505, Q12.97
	// (not the arithmetic Phred mean, 20).
	if got := MeanQScore([]byte{10, 30}, 50); math.Abs(got-(-10*math.Log10(0.0505))) > 1e-12 {
		t.Errorf("MeanQScore = %v", got)
	}
	if got := MeanQScore([]byte{60, 60}, 50); math.Abs(got-50) > 1e-9 {
		t.Errorf("scores above the cap: MeanQScore = %v, want 50", got)
	}
	if MeanQScore(nil, 50) != 0 {
		t.Error("empty read")
	}
}

func TestAppendSeconds(t *testing.T) {
	for d, want := range map[time.Duration]string{
		0: "0.000000", 1500 * time.Millisecond: "1.500000", 3*time.Hour + 7*time.Microsecond: "10800.000007",
	} {
		if got := string(appendSeconds(nil, d)); got != want {
			t.Errorf("appendSeconds(%v) = %q, want %q", d, got, want)
		}
	}
}
