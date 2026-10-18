package readlen

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Histogram maps a read length to its count (weight).
type Histogram map[int]float64

// Lengths returns the histogram's lengths in ascending order.
func (h Histogram) Lengths() []int {
	ls := make([]int, 0, len(h))
	for l := range h {
		ls = append(ls, l)
	}
	sort.Ints(ls)
	return ls
}

// LoadHistogram reads a read-length histogram file.
//
// Format: one "length count" pair per line, separated by a tab, a comma or
// spaces. length is a positive integer number of bases; count is a
// non-negative number (integer or decimal weight). Blank lines and lines
// starting with '#' are ignored, and so is a first non-comment line whose
// first field is not a number (a header such as "length<TAB>count").
// Repeated lengths are summed. The total count must be positive. A histogram
// can be derived from an ONT sequencing_summary.txt with
//
//	awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="sequence_length_template")c=i;next}
//	    {n[$c]++} END{for(l in n)print l"\t"n[l]}' sequencing_summary.txt
//
// Example:
//
//	# length	count
//	812	3
//	1450	11
//	23017	1
func LoadHistogram(path string) (Histogram, error) {
	f, err := os.Open(path) // #nosec G304 -- user-supplied histogram file is the point
	if err != nil {
		return nil, fmt.Errorf("%w: reading length histogram: %v", ErrInvalidSpec, err)
	}
	defer f.Close()

	h := Histogram{}
	var total float64
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 64*1024), 1024*1024)
	line, data := 0, 0
	for sc.Scan() {
		line++
		text := strings.TrimSpace(sc.Text())
		if text == "" || strings.HasPrefix(text, "#") {
			continue
		}
		data++
		fields := strings.FieldsFunc(text, func(r rune) bool { return r == '\t' || r == ',' || r == ' ' })
		if data == 1 && len(fields) > 0 && !numeric(fields[0]) {
			continue // header line
		}
		if len(fields) != 2 {
			return nil, fmt.Errorf("%w: %s:%d: want \"length<TAB>count\", got %q", ErrInvalidSpec, path, line, text)
		}
		length, err := strconv.Atoi(fields[0])
		if err != nil || length < 1 {
			return nil, fmt.Errorf("%w: %s:%d: length %q must be a positive integer", ErrInvalidSpec, path, line, fields[0])
		}
		count, err := strconv.ParseFloat(fields[1], 64)
		if err != nil || count < 0 || math.IsNaN(count) || math.IsInf(count, 0) {
			return nil, fmt.Errorf("%w: %s:%d: count %q must be a non-negative number", ErrInvalidSpec, path, line, fields[1])
		}
		h[length] += count
		total += count
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("%w: reading length histogram %s: %v", ErrInvalidSpec, path, err)
	}
	if !(total > 0) {
		return nil, fmt.Errorf("%w: length histogram %s has no counts", ErrInvalidSpec, path)
	}
	return h, nil
}

func numeric(s string) bool {
	_, err := strconv.ParseFloat(s, 64)
	return err == nil
}
