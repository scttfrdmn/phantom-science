// Phantom Science - Scientific Simulation Platform
// Copyright 2025-2026 Scott Friedman
// SPDX-License-Identifier: Apache-2.0

// Command phantom-science generates synthetic sequencing data (FASTQ and
// related formats) for testing pipelines. The subcommands and flags live in
// package cmd; run "phantom-science --help" for the list.
package main

import (
	"github.com/scttfrdmn/phantom-science/cmd"
)

func main() {
	cmd.Execute()
}
