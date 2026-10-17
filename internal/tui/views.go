package tui

import (
	"fmt"
	"path/filepath"
	"strings"
	"time"

	"github.com/charmbracelet/lipgloss"
	"github.com/scttfrdmn/phantom-science/internal/timing"
)

// renderHeader renders the main header
func (m Model) renderHeader() string {
	title := "Phantom Science Sequencer Simulator v0.3.0"
	return m.styles.Header.Render(title)
}

// renderTopPanels renders the sequencer selection, run status, and stats panels
func (m Model) renderTopPanels() string {
	sequencerPanel := m.renderSequencerPanel()
	statusPanel := m.renderStatusPanel()
	statsPanel := m.renderStatsPanel()

	return lipgloss.JoinHorizontal(
		lipgloss.Top,
		sequencerPanel,
		statusPanel,
		statsPanel,
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
	sequencers := []struct {
		name     string
		flowCell string
		icon     string
	}{
		{"NovaSeq 6000", "S2", "🚀"},
		{"NextSeq 2000", "P3", "⚡"},
		{"MiSeq", "v3", "🔬"},
	}

	for i, seq := range sequencers {
		prefix := "○"
		optionStyle := m.styles.SequencerOption
		
		if i == m.sequencerChoice {
			prefix = "●"
			optionStyle = m.styles.SequencerOptionActive
		}
		
		option := fmt.Sprintf("%s %s %s (%s)", prefix, seq.icon, seq.name, seq.flowCell)
		options = append(options, optionStyle.Render(option))
	}

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		strings.Join(options, "\n"),
	)

	return style.Width(25).Render(content)
}

// renderStatusPanel renders the run status panel
func (m Model) renderStatusPanel() string {
	style := m.styles.Panel
	
	title := m.styles.Subtitle.Render("Run Status")
	
	status := "Ready"
	progressPercent := 0.0
	eta := "--:--"
	
	if m.runInProgress {
		status = string(m.progress.Stage)
		progressPercent = m.progress.Percent
		if m.progress.ETA > 0 {
			eta = m.progress.ETA.Round(time.Second).String()
		}
	}

	statusLine := fmt.Sprintf("Status: %s", m.styles.StatusLabel(status))
	progressLine := fmt.Sprintf("Progress: %s %.1f%%", 
		m.styles.ProgressBarView(progressPercent, 10), progressPercent)
	etaLine := fmt.Sprintf("ETA: %s", m.styles.Value.Render(eta))

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

// renderStatsPanel renders the real-time statistics panel
func (m Model) renderStatsPanel() string {
	style := m.styles.Panel
	
	title := m.styles.Subtitle.Render("Real-time Stats")
	
	filesGenerated := "0"
	totalFiles := "0"
	totalSize := "0 MB"
	currentRate := "0 MB/s"
	currentFile := "None"
	
	if m.runInProgress {
		filesGenerated = fmt.Sprintf("%d", m.progress.FilesGenerated)
		totalFiles = fmt.Sprintf("%d", m.progress.TotalFiles)
		totalSize = fmt.Sprintf("%.1f MB", float64(m.progress.BytesGenerated)/(1024*1024))
		
		// Calculate current rate based on timing
		if m.simulator != nil {
			currentRate = fmt.Sprintf("%.1f MB/s", float64(m.simulator.Timing.BytesPerSecond)/(1024*1024))
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
	sizeLine := fmt.Sprintf("Size: %s", m.styles.Value.Render(totalSize))
	rateLine := fmt.Sprintf("Rate: %s", m.styles.Value.Render(currentRate))
	currentFileLine := fmt.Sprintf("Current: %s", m.styles.Value.Render(currentFile))

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

// renderSampleSheet renders the sample sheet panel
func (m Model) renderSampleSheet() string {
	style := m.styles.Panel
	if m.currentSection == SectionSamples {
		style = m.styles.ActivePanel
	}

	title := m.styles.Subtitle.Render("Sample Sheet")

	// Table header
	headers := []string{"Sample_ID", "Lane", "Index_1", "Index_2", "Project", "Description"}
	headerRow := m.styles.SampleHeader.Render(
		fmt.Sprintf("%-12s %-4s %-9s %-9s %-12s %-15s", headers[0], headers[1], headers[2], headers[3], headers[4], headers[5]))

	// Sample rows
	var rows []string
	for i, sample := range m.samples {
		rowStyle := m.styles.SampleRow
		if m.currentSection == SectionSamples && i == 0 { // Highlight first row when section is active
			rowStyle = m.styles.SampleRowActive
		}
		
		row := rowStyle.Render(
			fmt.Sprintf("%-12s %-4d %-9s %-9s %-12s %-15s", 
				sample.ID, sample.Lane, sample.Index1, sample.Index2, sample.Project, sample.Description))
		rows = append(rows, row)
	}

	// Action buttons
	buttons := lipgloss.JoinHorizontal(
		lipgloss.Left,
		m.styles.Button.Render("[Add Sample]"),
		m.styles.Button.Render("[Import CSV]"),
		m.styles.Button.Render("[Clear All]"),
	)

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		headerRow,
		strings.Join(rows, "\n"),
		"",
		buttons,
	)

	return style.Width(75).Render(content)
}

// renderBottomPanels renders the parameters, fingerprinting, and output settings panels
func (m Model) renderBottomPanels() string {
	parametersPanel := m.renderParametersPanel()
	fingerprintingPanel := m.renderFingerprintingPanel()
	outputPanel := m.renderOutputPanel()

	return lipgloss.JoinHorizontal(
		lipgloss.Top,
		parametersPanel,
		fingerprintingPanel,
		outputPanel,
	)
}

// renderParametersPanel renders the run parameters panel
func (m Model) renderParametersPanel() string {
	style := m.styles.Panel
	if m.currentSection == SectionParameters {
		style = m.styles.ActivePanel
	}

	title := m.styles.Subtitle.Render("Run Parameters")

	readLengthLine := fmt.Sprintf("%s %s bp",
		m.styles.Label.Render("Read Length:"),
		m.styles.Value.Render(fmt.Sprintf("[%d]", m.sequencerConfig.ReadLength)))

	pairedEndCheck := "✓"
	if !m.sequencerConfig.PairedEnd {
		pairedEndCheck = "✗"
	}
	pairedEndLine := fmt.Sprintf("%s [%s]",
		m.styles.Label.Render("Paired-End:"),
		m.styles.Value.Render(pairedEndCheck))

	chemistryLine := fmt.Sprintf("%s %s",
		m.styles.Label.Render("Chemistry:"),
		m.styles.Value.Render("[v1.5]"))

	phixLine := fmt.Sprintf("%s %s",
		m.styles.Label.Render("PhiX %:"),
		m.styles.Value.Render("[1.0]"))

	// HMM model — show basename so the path fits in the panel width
	hmmDisplay := "—"
	if m.hmmModelPath != "" {
		hmmDisplay = filepath.Base(m.hmmModelPath)
	}
	hmmLine := fmt.Sprintf("%s %s",
		m.styles.Label.Render("HMM Model:"),
		m.styles.Value.Render(fmt.Sprintf("[%s]", hmmDisplay)))

	// Population code (e.g. YRI, CEU, CHB)
	popDisplay := "—"
	if m.population != "" {
		popDisplay = m.population
	}
	populationLine := fmt.Sprintf("%s %s",
		m.styles.Label.Render("Population:"),
		m.styles.Value.Render(fmt.Sprintf("[%s]", popDisplay)))

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		readLengthLine,
		pairedEndLine,
		chemistryLine,
		phixLine,
		hmmLine,
		populationLine,
	)

	return style.Width(30).Render(content)
}

// renderFingerprintingPanel renders the enterprise fingerprinting panel
func (m Model) renderFingerprintingPanel() string {
	style := m.styles.Panel
	if m.currentSection == SectionFingerprinting {
		style = m.styles.ActivePanel
	}

	title := m.styles.Subtitle.Render("🔒 Enterprise Fingerprinting")

	// Fingerprinting enabled/disabled
	enabledCheck := "✗"
	if m.fingerprintingEnabled {
		enabledCheck = "✓"
	}
	enabledLine := fmt.Sprintf("%s [%s]", 
		m.styles.Label.Render("Enabled:"), 
		m.styles.Value.Render(enabledCheck))

	// Customer ID
	customerLine := fmt.Sprintf("%s %s", 
		m.styles.Label.Render("Customer ID:"), 
		m.styles.Value.Render(fmt.Sprintf("[%s]", m.customerID)))

	// Embedding strategy
	strategies := []string{"Conservative", "Balanced", "Aggressive"}
	strategyName := strategies[m.embeddingStrategy]
	strategyLine := fmt.Sprintf("%s %s", 
		m.styles.Label.Render("Strategy:"), 
		m.styles.Value.Render(fmt.Sprintf("[%s]", strategyName)))

	// Strategy details
	var strategyDetails string
	switch m.embeddingStrategy {
	case StrategyConservative:
		strategyDetails = m.styles.Subtitle.Render("1 bit/char, minimal risk")
	case StrategyBalanced:
		strategyDetails = m.styles.Subtitle.Render("2 bits/char, standard use")
	case StrategyAggressive:
		strategyDetails = m.styles.Subtitle.Render("3 bits/char, max capacity")
	}

	// Powered by stegano
	poweredBy := m.styles.Subtitle.Render("Powered by stegano v0.2.0")

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		enabledLine,
		customerLine,
		strategyLine,
		strategyDetails,
		"",
		poweredBy,
	)

	return style.Width(35).Render(content)
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
		
	compressionLine := fmt.Sprintf("%s [✓] gzip", 
		m.styles.Label.Render("Compression:"))
		
	realtimeLine := fmt.Sprintf("%s [✓] Stream files during run", 
		m.styles.Label.Render("Real-time:"))

	// Time compression options
	compressionOptions := []string{"[1x]", "[10x]", "[100x]", "[1000x]"}
	for i, option := range compressionOptions {
		if timing.TimeCompression((i+1)*10) == m.compression || (i == 0 && m.compression == 1) {
			compressionOptions[i] = m.styles.Value.Render(option)
		} else {
			compressionOptions[i] = m.styles.Value.Foreground(m.styles.Muted).Render(option)
		}
	}
	
	timeCompressionLine := fmt.Sprintf("%s %s", 
		m.styles.Label.Render("Time Compression:"), 
		strings.Join(compressionOptions, " "))

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		outputDirLine,
		compressionLine,
		realtimeLine,
		timeCompressionLine,
	)

	return style.Width(45).Render(content)
}

// renderPreview renders the file generation preview panel
func (m Model) renderPreview() string {
	style := m.styles.Panel
	
	title := m.styles.Subtitle.Render("File Generation Preview")

	// Calculate expected output
	totalFiles := m.sequencerConfig.TotalBatches
	if m.sequencerConfig.PairedEnd {
		totalFiles *= 2
	}
	
	// Estimate file size more accurately
	bytesPerRead := int64(50 + m.sequencerConfig.ReadLength*2 + 4) // Header + sequence + quality + newlines
	totalReads := int64(totalFiles * m.sequencerConfig.ReadsPerBatch)
	estimatedBytes := bytesPerRead * totalReads
	estimatedSizeGB := float64(estimatedBytes) / (1024 * 1024 * 1024)
	
	// Get timing estimate
	var duration string
	if m.simulator != nil {
		realHours := m.simulator.Timing.TotalRunTimeHours
		compressedMinutes := float64(realHours * 60) / float64(m.compression)
		if compressedMinutes < 60 {
			duration = fmt.Sprintf("~%d hours (real) / %.1f minutes (%dx compressed)", 
				realHours, compressedMinutes, m.compression)
		} else {
			duration = fmt.Sprintf("~%d hours (real) / %.1f hours (%dx compressed)", 
				realHours, compressedMinutes/60, m.compression)
		}
	} else {
		// Default estimates
		sequencerTiming := timing.GetSequencerTiming(m.sequencerConfig.SequencerModel, m.sequencerConfig.FlowCellType)
		realHours := sequencerTiming.TotalRunTimeHours
		compressedMinutes := float64(realHours * 60) / float64(m.compression)
		if compressedMinutes < 60 {
			duration = fmt.Sprintf("~%d hours (real) / %.1f minutes (%dx compressed)", 
				realHours, compressedMinutes, m.compression)
		} else {
			duration = fmt.Sprintf("~%d hours (real) / %.1f hours (%dx compressed)", 
				realHours, compressedMinutes/60, m.compression)
		}
	}

	// Generate sample filename
	sampleName := m.sequencerConfig.SampleName
	if len(m.samples) > 0 {
		sampleName = m.samples[0].ID
	}

	expectedLine := fmt.Sprintf("Expected Files: %s FASTQ files (~%.2f GB total)", 
		m.styles.Value.Render(fmt.Sprintf("%d", totalFiles)), estimatedSizeGB)
		
	namingLine := fmt.Sprintf("Naming: %s", 
		m.styles.Value.Render(fmt.Sprintf("%s_S1_L001_R1_001.fastq.gz", sampleName)))
		
	durationLine := fmt.Sprintf("Duration: %s", 
		m.styles.Value.Render(duration))
	
	// Output directory
	outputLine := fmt.Sprintf("Output: %s", 
		m.styles.Value.Render(m.sequencerConfig.OutputDir))

	// File generation rate
	var rateLine string
	if m.simulator != nil {
		filesPerHour := m.simulator.Timing.FilesPerHour
		rateLine = fmt.Sprintf("Rate: %s files/hour", 
			m.styles.Value.Render(fmt.Sprintf("%d", filesPerHour)))
	} else {
		sequencerTiming := timing.GetSequencerTiming(m.sequencerConfig.SequencerModel, m.sequencerConfig.FlowCellType)
		rateLine = fmt.Sprintf("Rate: %s files/hour", 
			m.styles.Value.Render(fmt.Sprintf("%d", sequencerTiming.FilesPerHour)))
	}

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

// renderActions renders the action buttons panel
func (m Model) renderActions() string {
	style := m.styles.Panel
	if m.currentSection == SectionActions {
		style = m.styles.ActivePanel
	}

	title := m.styles.Subtitle.Render("Actions")

	var buttons []string
	
	if m.runInProgress {
		buttons = append(buttons, m.styles.ButtonActive.Render("[Stop Run]"))
		buttons = append(buttons, m.styles.ButtonDisabled.Render("[Load Config]"))
		buttons = append(buttons, m.styles.ButtonDisabled.Render("[Save Config]"))
	} else {
		buttons = append(buttons, m.styles.ButtonActive.Render("[Start Run]"))
		buttons = append(buttons, m.styles.Button.Render("[Load Config]"))
		buttons = append(buttons, m.styles.Button.Render("[Save Config]"))
	}
	
	buttons = append(buttons, m.styles.Button.Render("[Help]"))
	buttons = append(buttons, m.styles.Button.Render("[Quit]"))

	buttonRow := lipgloss.JoinHorizontal(lipgloss.Left, buttons...)

	content := lipgloss.JoinVertical(
		lipgloss.Left,
		title,
		"",
		buttonRow,
	)

	return style.Width(75).Render(content)
}