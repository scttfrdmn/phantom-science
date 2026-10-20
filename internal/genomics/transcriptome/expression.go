package transcriptome

import (
	"fmt"
	"io"
	"math"
	"math/rand/v2"
	"strconv"
	"strings"
)

// Default expression-model parameters.
//
// Transcript abundance follows Zipf's law: the rank–frequency plot of
// expression is a power law with exponent close to −1 across cell types and
// tissues (Furusawa & Kaneko, "Zipf's law in gene expression", Phys Rev Lett
// 90:088102, 2003; Ramsköld et al., PLoS Comput Biol 5:e1000598, 2009, for
// RNA-seq, whose expression distributions are also commonly described as
// log-normal). With α = 1 and no scatter, the top k ≈ e^{(ln n − γ)/2} ≈ 375
// of GENCODE's ≈250,000 transcripts carry half of all molecules
// (H(k) = H(n)/2); with the default scatter the measured figure is ≈250.
// That is the "few hundred transcripts carry half the reads" of real brain
// cDNA. Each transcript's level is scattered log-normally around its rank's
// power-law value (log-normal given rank). The scatter σ = 0.5 (a factor of
// ≈1.6) is a modelling choice.
const (
	DefaultZipfAlpha       = 1.0
	DefaultExpressionSigma = 0.5
)

// ZipfModel draws transcript abundances: transcripts get a random rank r
// (a permutation of 1..n) and
//
//	log TPM_i = c − Alpha·log r_i + Sigma·Z_i,  Z_i ~ N(0, 1),
//
// with c normalising the TPMs to sum to 10^6.
type ZipfModel struct {
	Alpha float64
	Sigma float64
}

// Validate checks the parameters.
func (m ZipfModel) Validate() error {
	if !(m.Alpha >= 0) || math.IsInf(m.Alpha, 0) {
		return fmt.Errorf("zipf alpha must be finite and >= 0, got %g", m.Alpha)
	}
	if !(m.Sigma >= 0) || math.IsInf(m.Sigma, 0) {
		return fmt.Errorf("expression sigma must be finite and >= 0, got %g", m.Sigma)
	}
	return nil
}

// TPM draws n abundances from r. It returns the TPMs (summing to 10^6) and
// the 1-based rank assigned to each transcript.
func (m ZipfModel) TPM(r *rand.Rand, n int) (tpm []float64, rank []int) {
	perm := r.Perm(n)
	rank = make([]int, n)
	logw := make([]float64, n)
	mx := math.Inf(-1)
	for i := range n {
		rank[i] = perm[i] + 1
		logw[i] = float64(-m.Alpha*math.Log(float64(rank[i]))) + float64(m.Sigma*r.NormFloat64())
		mx = max(mx, logw[i])
	}
	tpm = make([]float64, n)
	sum := 0.0
	for i, lw := range logw {
		tpm[i] = math.Exp(lw - mx)
		sum += tpm[i]
	}
	for i := range tpm {
		tpm[i] = float64(tpm[i] * (1e6 / sum))
	}
	return tpm, rank
}

// ProfileStats reports how an expression profile matched the table.
type ProfileStats struct {
	Rows      int // data rows read
	Matched   int // rows matched to a table transcript
	Unmatched int // rows whose ID is not in the table
	Column    string
}

// LoadProfile reads an expression profile (plain or gzip) and returns one
// TPM per table transcript (summing to 10^6; transcripts absent from the
// profile get 0). See ParseProfile for the formats.
func LoadProfile(path string, t *Table) ([]float64, ProfileStats, error) {
	rc, err := openInput(path)
	if err != nil {
		return nil, ProfileStats{}, err
	}
	defer rc.Close() // #nosec G307 -- read-only
	tpm, st, err := ParseProfile(rc, t)
	if err != nil {
		return nil, st, fmt.Errorf("%s: %w", path, err)
	}
	return tpm, st, nil
}

// ParseProfile reads a tab-separated expression table. Accepted formats:
//
//   - two columns "transcript_id<TAB>tpm", with or without a header line;
//   - a table with a header naming a TPM column ("TPM"/"tpm", e.g. salmon
//     quant.sf, which has Name/Length/EffectiveLength/TPM/NumReads);
//   - Bambu counts (counts_transcript.txt / transcript_counts.tsv): header
//     "TXNAME<TAB>GENEID<TAB><sample>…"; the sample columns are summed.
//
// The ID column is the header column named transcript_id, TXNAME, Name,
// target_id or tx_id, else the first. Lines starting with '#' are ignored.
// Values are relative abundances: they are rescaled to sum to 10^6, so
// counts work as well as TPM. For full-length cDNA one read is one molecule,
// so read counts are proportional to TPM without length correction.
//
// IDs match exactly, or, failing that, with the version suffix (".N")
// removed from both sides (ENST00000456328.2 matches ENST00000456328). Rows
// for the same transcript are summed.
func ParseProfile(r io.Reader, t *Table) ([]float64, ProfileStats, error) {
	var st ProfileStats
	unversioned := map[string]int{}
	for i, tx := range t.Transcripts {
		k := stripVersion(tx.ID)
		if _, dup := unversioned[k]; dup {
			unversioned[k] = -1 // ambiguous
		} else {
			unversioned[k] = i
		}
	}
	tpm := make([]float64, len(t.Transcripts))
	lr := newLineReader(r)
	idCol, valCols := 0, []int{1}
	first := true
	for {
		line, err := lr.next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, st, err
		}
		s := strings.TrimSpace(string(line))
		if s == "" || s[0] == '#' {
			continue
		}
		f := strings.Split(s, "\t")
		if first {
			first = false
			if hdr, ok := profileHeader(f); ok {
				idCol, valCols, st.Column = hdr.id, hdr.vals, hdr.name
				continue
			}
			st.Column = "column 2"
		}
		need := idCol
		for _, c := range valCols {
			need = max(need, c)
		}
		if len(f) <= need {
			return nil, st, fmt.Errorf("line %d: want at least %d columns, got %d", lr.line, need+1, len(f))
		}
		v := 0.0
		for _, c := range valCols {
			x, err := strconv.ParseFloat(strings.TrimSpace(f[c]), 64)
			if err != nil || x < 0 || math.IsInf(x, 0) || math.IsNaN(x) {
				return nil, st, fmt.Errorf("line %d: bad abundance %q", lr.line, f[c])
			}
			v += x
		}
		st.Rows++
		id := strings.TrimSpace(f[idCol])
		i := t.Index(id)
		if i < 0 {
			if j, ok := unversioned[stripVersion(id)]; ok && j >= 0 {
				i = j
			}
		}
		if i < 0 {
			st.Unmatched++
			continue
		}
		st.Matched++
		tpm[i] += v
	}
	sum := 0.0
	for _, v := range tpm {
		sum += v
	}
	if st.Matched == 0 {
		return nil, st, fmt.Errorf("no profile ID matches a transcript (%d rows)", st.Rows)
	}
	if !(sum > 0) {
		return nil, st, fmt.Errorf("matched transcripts have zero total abundance")
	}
	for i := range tpm {
		tpm[i] = float64(tpm[i] * (1e6 / sum))
	}
	return tpm, st, nil
}

type headerCols struct {
	id   int
	vals []int
	name string
}

// profileHeader recognises a header row: one whose value column(s) are not
// numeric.
func profileHeader(f []string) (headerCols, bool) {
	if len(f) < 2 {
		return headerCols{}, false
	}
	if _, err := strconv.ParseFloat(strings.TrimSpace(f[1]), 64); err == nil {
		return headerCols{}, false
	}
	h := headerCols{id: 0}
	lower := make([]string, len(f))
	for i, c := range f {
		lower[i] = strings.ToLower(strings.TrimSpace(c))
	}
	for i, c := range lower {
		switch c {
		case "transcript_id", "txname", "name", "target_id", "tx_id":
			h.id = i
		}
	}
	for i, c := range lower {
		if c == "tpm" {
			h.vals, h.name = []int{i}, f[i]
			return h, true
		}
	}
	if len(f) >= 3 && lower[0] == "txname" && lower[1] == "geneid" {
		for i := 2; i < len(f); i++ {
			h.vals = append(h.vals, i)
		}
		h.name = "sum of Bambu sample counts"
		return h, true
	}
	for i := range f {
		if i != h.id {
			h.vals, h.name = []int{i}, f[i]
			return h, true
		}
	}
	return headerCols{}, false
}

// stripVersion removes a trailing ".<digits>" version from an Ensembl ID.
func stripVersion(id string) string {
	i := strings.LastIndexByte(id, '.')
	if i <= 0 || i == len(id)-1 {
		return id
	}
	for _, c := range id[i+1:] {
		if c < '0' || c > '9' {
			return id
		}
	}
	return id[:i]
}
