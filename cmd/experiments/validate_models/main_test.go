//go:build !release
// +build !release

package main

import (
	"testing"
)

func TestFindMaxHomopolymerRun_Empty(t *testing.T) {
	got := findMaxHomopolymerRun([]byte{})
	if got != 0 {
		t.Errorf("got %v, want 0", got)
	}
}

func TestFindMaxHomopolymerRun_SingleBase(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("A"))
	if got != 1 {
		t.Errorf("got %v, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_NoRun(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("ACGT"))
	if got != 1 {
		t.Errorf("got %v, want 1 for alternating sequence", got)
	}
}

func TestFindMaxHomopolymerRun_RunOf4(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("ACGTAAAACGT"))
	if got != 4 {
		t.Errorf("got %v, want 4", got)
	}
}

func TestFindMaxHomopolymerRun_RunAtEnd(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("ACGTTTTT"))
	if got != 5 {
		t.Errorf("got %v, want 5 for run at end", got)
	}
}

func TestFindMaxHomopolymerRun_RunAtStart(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("CCCACGT"))
	if got != 3 {
		t.Errorf("got %v, want 3 for run at start", got)
	}
}

func TestFindMaxHomopolymerRun_AllSame(t *testing.T) {
	seq := []byte("AAAAAAAAAA")
	got := findMaxHomopolymerRun(seq)
	if got != len(seq) {
		t.Errorf("got %v, want %d for all-same sequence", got, len(seq))
	}
}

func TestFindMaxHomopolymerRun_TwoRuns_TakesMax(t *testing.T) {
	// GGG (3) then TTTTT (5) -> max is 5
	got := findMaxHomopolymerRun([]byte("ACGGGACGTTTTTACG"))
	if got != 5 {
		t.Errorf("got %v, want 5", got)
	}
}

func TestFindMaxHomopolymerRun_TwoBases(t *testing.T) {
	// Two-character sequence with same bases
	got := findMaxHomopolymerRun([]byte("AA"))
	if got != 2 {
		t.Errorf("got %v, want 2", got)
	}
}

func TestFindMaxHomopolymerRun_TwoDifferentBases(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("AT"))
	if got != 1 {
		t.Errorf("got %v, want 1", got)
	}
}

func TestFindMaxHomopolymerRun_MultipleSameRunLength(t *testing.T) {
	// AAA and CCC are both length 3; result should be 3
	got := findMaxHomopolymerRun([]byte("AAATCCC"))
	if got != 3 {
		t.Errorf("got %v, want 3", got)
	}
}

func TestFindMaxHomopolymerRun_LongSequenceNoRun(t *testing.T) {
	got := findMaxHomopolymerRun([]byte("ACGTACGTACGTACGT"))
	if got != 1 {
		t.Errorf("got %v, want 1 for fully alternating sequence", got)
	}
}

func TestFindMaxHomopolymerRun_RunOf10InMiddle(t *testing.T) {
	// 10 G's in the middle
	seq := []byte("ACGT" + "GGGGGGGGGG" + "ACGT")
	got := findMaxHomopolymerRun(seq)
	if got != 10 {
		t.Errorf("got %v, want 10", got)
	}
}

func TestFindMaxHomopolymerRun_TableDriven(t *testing.T) {
	cases := []struct {
		name  string
		seq   string
		want  int
	}{
		{"empty", "", 0},
		{"single A", "A", 1},
		{"single G", "G", 1},
		{"two same", "TT", 2},
		{"two diff", "AC", 1},
		{"all C run of 5", "CCCCC", 5},
		{"run at start 2", "GGAC", 2},
		{"run at end 2", "ACGG", 2},
		{"long alternating", "ATATATATAT", 1},
		{"mixed runs pick largest", "AACCCTTTT", 4},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := findMaxHomopolymerRun([]byte(tc.seq))
			if got != tc.want {
				t.Errorf("findMaxHomopolymerRun(%q) = %d, want %d", tc.seq, got, tc.want)
			}
		})
	}
}
