package transcriptome

import (
	"fmt"
	"sort"
	"strings"
)

// Kit holds the library-preparation sequences that flank each cDNA molecule.
//
// A sense ("+") read is
//
//	Adapter5 · SSP · transcript (3′ suffix) · poly(A) · revcomp(VNP) · Adapter3
//
// and a "−" read is the reverse complement of the whole construct. This is
// the read structure Pychopper classifies (SSP at the 5′ end and the
// reverse-complemented VNP at the 3′ end of forward reads).
//
// SSP and VNP are the strand-switching and RT (oligo-dT VN) primer
// sequences of Pychopper's kit primer files
// (https://github.com/epi2me-labs/pychopper, pychopper/primer_data:
// cDNA_SSP_VNP.fas for PCS109, PCS110_primers.fas, PCS111_primers.fas),
// with the oligo-dT tail of VNP removed: the poly(A) drawn for each molecule
// stands for the tail and the oligo-dT that pairs with it. Pychopper has no
// separate PCS114 file and its kit choices are PCS109/PCS110/PCS111/LSK114;
// PCS114 is the Kit 14 release of PCS111 and is assumed here to use the
// PCS111 primers (run Pychopper with -k PCS111). This is an assumption, not
// a published statement: override the sequences if your kit differs.
//
// Adapter5 and Adapter3 are empty by default. Dorado trims the sequencing
// adapter by default, so basecalled reads usually start near the primer;
// set them to model untrimmed adapter sequence.
type Kit struct {
	Name     string
	SSP      string
	VNP      string
	Adapter5 string
	Adapter3 string
}

// Published primer sequences (Pychopper primer_data; oligo-dT removed).
const (
	sspPCS109 = "TTTCTGTTGGTGCTGATATTGCTGGG"
	vnpPCS109 = "ACTTGCCTGTCGCTCTATCTTC"
	sspPCS111 = "TTTCTGTTGGTGCTGATATTGCTTT"
	vnpPCS110 = "ACTTGCCTGTCGCTCTATCTTCGAGGAGAGTCCGCCGCCCGCAAG"
	vnpPCS111 = "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG"
)

// Kits are the built-in cDNA-PCR kits.
var Kits = map[string]Kit{
	"PCS109": {Name: "PCS109", SSP: sspPCS109, VNP: vnpPCS109},
	"PCS110": {Name: "PCS110", SSP: sspPCS111, VNP: vnpPCS110},
	"PCS111": {Name: "PCS111", SSP: sspPCS111, VNP: vnpPCS111},
	"PCS114": {Name: "PCS114", SSP: sspPCS111, VNP: vnpPCS111},
}

// DefaultKit is the current ONT cDNA-PCR kit.
const DefaultKit = "PCS114"

// KitNames returns the built-in kit names, sorted.
func KitNames() []string {
	var n []string
	for k := range Kits {
		n = append(n, k)
	}
	sort.Strings(n)
	return n
}

// LookupKit returns the built-in kit name (case-insensitive).
func LookupKit(name string) (Kit, error) {
	k, ok := Kits[strings.ToUpper(name)]
	if !ok {
		return Kit{}, fmt.Errorf("unknown cDNA kit %q: want one of %s", name, strings.Join(KitNames(), ", "))
	}
	return k, nil
}

func (k Kit) validate() error {
	for _, f := range []struct{ name, s string }{
		{"SSP", k.SSP}, {"VNP", k.VNP}, {"Adapter5", k.Adapter5}, {"Adapter3", k.Adapter3},
	} {
		if i := strings.IndexFunc(f.s, func(c rune) bool { return !strings.ContainsRune("ACGT", c) }); i >= 0 {
			return fmt.Errorf("kit %s %s contains %q: only ACGT allowed", k.Name, f.name, f.s[i])
		}
	}
	return nil
}
