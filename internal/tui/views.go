package tui

import (
	"fmt"
	"path/filepath"
	"strings"
	"time"

	"github.com/charmbracelet/lipgloss"
	"github.com/scttfrdmn/phantom-science/internal/timing"
)

// orDash returns s, or "—" when s is empty.
func orDash(s string) string {
	if s == "" {
		return "—"
	}
	return s
}

// yesNo renders a boolean setting.
func yesNo(b bool) string {
	if b {
		return "✓"
	}
	return "✗"
}

// runStarted reports whether a run has been started in this session.
func (m Model) runStarted() bool { return m.run != nil }

// renderHeader renders the main header
func (m Model) renderHeader() string {
	title := "Phantom Science Sequencer Simulator"
	if m.version != "" {
		title += " v" + strings.TrimPrefix(m.version, "v")
	}
	return m.styles.Header.Render(title)
}

// renderError renders the last run error.
func (m Model) renderError() string {
	msg := lipgloss.NewStyle().Foreground(m.styles.Error).Bold(true).Render("Error: ") + m.err.Error()
	return m.styles.Panel.Width(75).Render(msg)
}

// renderTopPanels renders the sequencer selection, run status, and stats panels
func (m Model) renderTopPanels() string {
	return lipgloss.JoinHorizontal(
		lipgloss.Top,
		m.renderSequencerPanel(),
		m.renderStatusPanel(),
		m.renderStatsPanel(),
	)
}

// renderSequencerPanel renders the sequencer selection panel
func (m Model) renderSequencerPanel() string {
	style := m.styles.Panel
	if m.currentSection == SectionSequencer {
		style = m.styles.ActivePanel
	}

	title := m.styles.Subtitle.Render("Sequencer Selection")

	var options []string
	if m.sequencerChoice == ChoiceConfig {
		// The configured sequencer is not a preset; it is used as configured.
		option := fmt.Sprintf("● %s (%s)",
			orDash(m.sequencerConfig.SequencerModel), orDash(m.sequencerConfig.FlowCellType))
		options = append(options, m.styles.SequencerOptionActive.Render(option))
	}
	for i, seq := range sequencerPresets {
		prefix := "○"
		optionStyle := m.styles.SequencerOption
		flowCell := seq.flowCell

		if i == m.sequencerChoice {
			prefix = "●"
			optionStyle = m.styles.SequencerOptionActive
			if !m.sequencerChanged {
				flowCell = orDash(m.sequencerConfig.FlowCellType) // as configured
			}
		}

		option := fmt.Sprintf("%s %s %s (%s)", prefix, seq.icon, seq.name, flowCell)
		options = append(options, optionStyle.Render(option))
	}

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		strings.Join(options, "\n"),
	)

	return style.Width(30).Render(content)
}

// renderStatusPanel renders the run status panel
func (m Model) renderStatusPanel() string {
	style := m.styles.Panel

	title := m.styles.Subtitle.Render("Run Status")

	status := "Ready"
	progressPercent := 0.0
	eta := "--:--"

	if m.runStarted() {
		status = string(m.progress.Stage)
		progressPercent = m.progress.Percent
		if m.runInProgress && m.progress.ETA > 0 {
			eta = m.progress.ETA.Round(time.Second).String()
		}
	}

	statusLine := fmt.Sprintf("Status: %s", m.styles.StatusLabel(status))
	progressLine := fmt.Sprintf("Stage: %s %.1f%%",
		m.styles.ProgressBarView(progressPercent, 10), progressPercent)
	etaLine := fmt.Sprintf("Stage ETA: %s", m.styles.Value.Render(eta))

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		statusLine,
		progressLine,
		etaLine,
	)

	return style.Width(25).Render(content)
}

// renderStatsPanel renders the files written by the current or last run.
func (m Model) renderStatsPanel() string {
	style := m.styles.Panel

	title := m.styles.Subtitle.Render("Run Output")

	filesGenerated := "0"
	totalFiles := "0"
	totalSize := "0.0 MB"
	currentRate := "—"
	currentFile := "None"

	if m.runStarted() {
		filesGenerated = fmt.Sprintf("%d", m.progress.FilesGenerated)
		totalFiles = fmt.Sprintf("%d", m.progress.TotalFiles)
		totalSize = fmt.Sprintf("%.1f MB", float64(m.progress.BytesGenerated)/(1024*1024))

		// Measured: bytes written so far over wall-clock time since start.
		if _, elapsed, _ := m.run.snapshot(); m.runInProgress && elapsed > 0 {
			currentRate = fmt.Sprintf("%.2f MB/s",
				float64(m.progress.BytesGenerated)/(1024*1024)/elapsed.Seconds())
		}

		if m.progress.CurrentFile != "" {
			// Show just the filename, not full path
			currentFile = filepath.Base(m.progress.CurrentFile)
			if len(currentFile) > 20 {
				currentFile = currentFile[:17] + "..."
			}
		}
	}

	filesLine := fmt.Sprintf("Files: %s/%s",
		m.styles.Value.Render(filesGenerated),
		m.styles.Value.Foreground(m.styles.Muted).Render(totalFiles))
	sizeLine := fmt.Sprintf("Written: %s", m.styles.Value.Render(totalSize))
	rateLine := fmt.Sprintf("Rate: %s", m.styles.Value.Render(currentRate))
	currentFileLine := fmt.Sprintf("Last: %s", m.styles.Value.Render(currentFile))

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		filesLine,
		sizeLine,
		rateLine,
		currentFileLine,
	)

	return style.Width(25).Render(content)
}

// renderSampleSheet renders the sample the run generates (sample_name and
// lane from the configuration). It is display-only.
func (m Model) renderSampleSheet() string {
	style := m.styles.Panel
	if m.currentSection == SectionSamples {
		style = m.styles.ActivePanel
	}

	title := m.styles.Subtitle.Render("Sample Sheet (from config)")

	headers := []string{"Sample_ID", "Lane", "Index_1", "Index_2", "Project", "Description"}
	headerRow := m.styles.SampleHeader.Render(
		fmt.Sprintf("%-12s %-4s %-9s %-9s %-12s %-15s", headers[0], headers[1], headers[2], headers[3], headers[4], headers[5]))

	var rows []string
	for _, sample := range m.samples {
		row := m.styles.SampleRow.Render(
			fmt.Sprintf("%-12s %-4d %-9s %-9s %-12s %-15s",
				orDash(sample.ID), sample.Lane, orDash(sample.Index1), orDash(sample.Index2),
				orDash(sample.Project), orDash(sample.Description)))
		rows = append(rows, row)
	}

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		headerRow,
		strings.Join(rows, "\n"),
	)

	return style.Width(75).Render(content)
}

// renderBottomPanels renders the parameters, fingerprinting, and output settings panels
func (m Model) renderBottomPanels() string {
	return lipgloss.JoinHorizontal(
		lipgloss.Top,
		m.renderParametersPanel(),
		m.renderFingerprintingPanel(),
		m.renderOutputPanel(),
	)
}

// renderParametersPanel renders the run parameters from the configuration.
func (m Model) renderParametersPanel() string {
	style := m.styles.Panel
	if m.currentSection == SectionParameters {
		style = m.styles.ActivePanel
	}
	c := m.sequencerConfig

	line := func(label, value string) string {
		return fmt.Sprintf("%s %s", m.styles.Label.Render(label), m.styles.Value.Render("["+value+"]"))
	}

	platform := c.Platform
	if platform == "" {
		platform = "illumina"
	}

	seed := "random"
	if c.Seed != 0 {
		seed = fmt.Sprintf("%d", c.Seed)
	}

	// HMM model — show basename so the path fits in the panel width
	hmmDisplay := "—"
	if c.HMMModelPath != "" {
		hmmDisplay = filepath.Base(c.HMMModelPath)
	}

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		m.styles.Subtitle.Render("Run Parameters"),
		"",
		line("Platform:", platform),
		line("Read Length:", fmt.Sprintf("%d bp", c.ReadLength)),
		line("Paired-End:", yesNo(c.PairedEnd)),
		line("Reads/Batch:", fmt.Sprintf("%d", c.ReadsPerBatch)),
		line("Batches:", fmt.Sprintf("%d", c.TotalBatches)),
		line("Seed:", seed),
		line("HMM Model:", hmmDisplay),
		line("Population:", orDash(c.Population)),
	)

	return style.Width(34).Render(content)
}

// renderFingerprintingPanel renders the fingerprinting settings the run uses.
func (m Model) renderFingerprintingPanel() string {
	style := m.styles.Panel
	if m.currentSection == SectionFingerprinting {
		style = m.styles.ActivePanel
	}

	layers := "—"
	if len(m.fingerprintLayers) > 0 {
		layers = strings.Join(m.fingerprintLayers, ",")
	}

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		m.styles.Subtitle.Render("🔒 Fingerprinting"),
		"",
		fmt.Sprintf("%s [%s]", m.styles.Label.Render("Enabled:"), m.styles.Value.Render(yesNo(m.fingerprintingEnabled))),
		fmt.Sprintf("%s %s", m.styles.Label.Render("Customer ID:"), m.styles.Value.Render("["+orDash(m.customerID)+"]")),
		fmt.Sprintf("%s %s", m.styles.Label.Render("Layers:"), m.styles.Value.Render("["+layers+"]")),
	)

	return style.Width(35).Render(content)
}

// outputCompression is the configured output compression ("none" when unset).
func (m Model) outputCompression() string {
	if m.sequencerConfig.Output.Compression == "" {
		return "none"
	}
	return m.sequencerConfig.Output.Compression
}

// renderOutputPanel renders the output settings panel
func (m Model) renderOutputPanel() string {
	style := m.styles.Panel
	if m.currentSection == SectionOutput {
		style = m.styles.ActivePanel
	}

	title := m.styles.Subtitle.Render("Output Settings")

	outputDirLine := fmt.Sprintf("%s %s",
		m.styles.Label.Render("Output Dir:"),
		m.styles.Value.Render(fmt.Sprintf("[%s]", m.sequencerConfig.OutputDir)))

	compressionLine := fmt.Sprintf("%s %s",
		m.styles.Label.Render("Compression:"),
		m.styles.Value.Render("["+m.outputCompression()+"]"))

	// Time compression options; the selected one is bracketed (←/→ changes it)
	options := make([]string, len(compressionOptions))
	for i, c := range compressionOptions {
		if c == m.compression {
			options[i] = m.styles.Value.Render(fmt.Sprintf("[%dx]", c))
		} else {
			options[i] = m.styles.Value.Foreground(m.styles.Muted).Render(fmt.Sprintf(" %dx ", c))
		}
	}

	timeCompressionLine := fmt.Sprintf("%s %s",
		m.styles.Label.Render("Time Speedup:"),
		strings.Join(options, " "))

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		outputDirLine,
		compressionLine,
		timeCompressionLine,
	)

	return style.Width(45).Render(content)
}

// renderPreview renders the output a run started now would produce.
func (m Model) renderPreview() string {
	style := m.styles.Panel

	title := m.styles.Subtitle.Render("File Generation Preview")

	c := m.updateConfigFromUI().sequencerConfig // what a run started now would use
	sim := timing.NewRunSimulator(c, m.compression)
	totalFiles, totalBytes := sim.EstimateTotalOutput()
	estimatedSize := fmt.Sprintf("%.1f MB", float64(totalBytes)/(1024*1024))
	if totalBytes >= 1024*1024*1024 {
		estimatedSize = fmt.Sprintf("%.2f GB", float64(totalBytes)/(1024*1024*1024))
	}

	realHours := sim.Timing.TotalRunTimeHours
	compressedMinutes := float64(realHours*60) / float64(m.compression)
	var duration string
	if compressedMinutes < 60 {
		duration = fmt.Sprintf("~%d hours (real) / %.1f minutes (%dx compressed)",
			realHours, compressedMinutes, m.compression)
	} else {
		duration = fmt.Sprintf("~%d hours (real) / %.1f hours (%dx compressed)",
			realHours, compressedMinutes/60, m.compression)
	}

	// The simulator names files <sample>_S1_L<lane>_R<mate>_<batch>.fastq; the
	// engine appends .gz for gzip output.
	firstFile := fmt.Sprintf("%s_S1_L%03d_R1_%03d.fastq", c.SampleName, c.Lane, 1)
	if strings.EqualFold(c.Output.Compression, "gzip") {
		firstFile += ".gz"
	}

	expectedLine := fmt.Sprintf("Expected Files: %s FASTQ files (~%s uncompressed)",
		m.styles.Value.Render(fmt.Sprintf("%d", totalFiles)), estimatedSize)
	namingLine := fmt.Sprintf("Naming: %s", m.styles.Value.Render(firstFile))
	outputLine := fmt.Sprintf("Output: %s", m.styles.Value.Render(c.OutputDir))
	rateLine := fmt.Sprintf("Instrument: %s",
		m.styles.Value.Render(fmt.Sprintf("%s (%s)", orDash(c.SequencerModel), orDash(c.FlowCellType))))
	durationLine := fmt.Sprintf("Duration: %s", m.styles.Value.Render(duration))

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		expectedLine,
		namingLine,
		outputLine,
		rateLine,
		durationLine,
	)

	return style.Width(75).Render(content)
}

// renderActions renders the start/stop button and the key bindings.
func (m Model) renderActions() string {
	style := m.styles.Panel
	if m.currentSection == SectionActions {
		style = m.styles.ActivePanel
	}

	title := m.styles.Subtitle.Render("Actions")

	button := m.styles.ButtonActive.Render("[Start Run]")
	if m.runInProgress {
		button = m.styles.ButtonActive.Render("[Stop Run]")
	}

	keys := m.styles.Value.Foreground(m.styles.Muted).Render(
		"Tab/Shift+Tab: panel  q: quit\n" +
			"Sequencer panel ↑/↓: sequencer  Output panel ←/→: time speedup\n" +
			"Actions panel Enter: start/stop run")

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		button,
		"",
		keys,
	)

	return style.Width(75).Render(content)
}
