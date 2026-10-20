package tui

import (
	"github.com/charmbracelet/lipgloss"
)

// Styles holds all the styling for the TUI
type Styles struct {
	// Main layout
	MainBorder  lipgloss.Style
	Header      lipgloss.Style
	Panel       lipgloss.Style
	ActivePanel lipgloss.Style

	// Sequencer selection
	SequencerOption       lipgloss.Style
	SequencerOptionActive lipgloss.Style

	// Sample sheet
	SampleTable     lipgloss.Style
	SampleHeader    lipgloss.Style
	SampleRow       lipgloss.Style
	SampleRowActive lipgloss.Style

	// Progress and status
	ProgressBar    lipgloss.Style
	ProgressFilled lipgloss.Style
	ProgressEmpty  lipgloss.Style
	StatusReady    lipgloss.Style
	StatusRunning  lipgloss.Style
	StatusComplete lipgloss.Style

	// Actions
	Button         lipgloss.Style
	ButtonActive   lipgloss.Style
	ButtonDisabled lipgloss.Style

	// Text elements
	Title            lipgloss.Style
	Subtitle         lipgloss.Style
	Label            lipgloss.Style
	Value            lipgloss.Style
	InputField       lipgloss.Style
	InputFieldActive lipgloss.Style

	// Colors
	Primary   lipgloss.Color
	Secondary lipgloss.Color
	Success   lipgloss.Color
	Warning   lipgloss.Color
	Error     lipgloss.Color
	Muted     lipgloss.Color
}

// NewStyles creates a new Styles instance with default styling
func NewStyles() Styles {
	// Color palette
	primary := lipgloss.Color("#7C3AED")    // Purple
	secondary := lipgloss.Color("#06B6D4")  // Cyan
	success := lipgloss.Color("#10B981")    // Green
	warning := lipgloss.Color("#F59E0B")    // Yellow
	errorColor := lipgloss.Color("#EF4444") // Red
	muted := lipgloss.Color("#6B7280")      // Gray

	return Styles{
		Primary:   primary,
		Secondary: secondary,
		Success:   success,
		Warning:   warning,
		Error:     errorColor,
		Muted:     muted,

		// Main layout
		MainBorder: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(primary).
			Padding(1).
			Margin(1),

		Header: lipgloss.NewStyle().
			Foreground(primary).
			Bold(true).
			Align(lipgloss.Center).
			Padding(0, 1),

		Panel: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(muted).
			Padding(1).
			Margin(0, 1),

		ActivePanel: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(primary).
			Padding(1).
			Margin(0, 1),

		// Sequencer selection
		SequencerOption: lipgloss.NewStyle().
			Padding(0, 2).
			Foreground(muted),

		SequencerOptionActive: lipgloss.NewStyle().
			Padding(0, 2).
			Foreground(primary).
			Bold(true),

		// Sample sheet
		SampleTable: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(muted).
			Width(70),

		SampleHeader: lipgloss.NewStyle().
			Foreground(primary).
			Bold(true).
			Padding(0, 1).
			BorderBottom(true).
			BorderStyle(lipgloss.NormalBorder()).
			BorderForeground(muted),

		SampleRow: lipgloss.NewStyle().
			Padding(0, 1).
			Foreground(lipgloss.Color("#FFFFFF")),

		SampleRowActive: lipgloss.NewStyle().
			Padding(0, 1).
			Foreground(primary).
			Background(lipgloss.Color("#1F2937")),

		// Progress and status
		ProgressBar: lipgloss.NewStyle().
			Width(20).
			Height(1),

		ProgressFilled: lipgloss.NewStyle().
			Background(success),

		ProgressEmpty: lipgloss.NewStyle().
			Background(muted),

		StatusReady: lipgloss.NewStyle().
			Foreground(success).
			Bold(true),

		StatusRunning: lipgloss.NewStyle().
			Foreground(warning).
			Bold(true),

		StatusComplete: lipgloss.NewStyle().
			Foreground(success).
			Bold(true),

		// Actions
		Button: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(muted).
			Padding(0, 2).
			Margin(0, 1).
			Foreground(lipgloss.Color("#FFFFFF")),

		ButtonActive: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(primary).
			Background(primary).
			Padding(0, 2).
			Margin(0, 1).
			Foreground(lipgloss.Color("#FFFFFF")).
			Bold(true),

		ButtonDisabled: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(muted).
			Padding(0, 2).
			Margin(0, 1).
			Foreground(muted),

		// Text elements
		Title: lipgloss.NewStyle().
			Foreground(primary).
			Bold(true).
			Padding(0, 0, 1, 0),

		Subtitle: lipgloss.NewStyle().
			Foreground(secondary).
			Bold(true),

		Label: lipgloss.NewStyle().
			Foreground(muted).
			Width(15),

		Value: lipgloss.NewStyle().
			Foreground(lipgloss.Color("#FFFFFF")).
			Bold(true),

		InputField: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(muted).
			Padding(0, 1).
			Width(20),

		InputFieldActive: lipgloss.NewStyle().
			Border(lipgloss.RoundedBorder()).
			BorderForeground(primary).
			Padding(0, 1).
			Width(20).
			Background(lipgloss.Color("#1F2937")),
	}
}

// ProgressBarView renders a progress bar
func (s Styles) ProgressBarView(percent float64, width int) string {
	if width <= 0 {
		width = 20
	}

	filled := int(percent * float64(width) / 100)
	if filled > width {
		filled = width
	}

	var bar string
	for i := 0; i < width; i++ {
		if i < filled {
			bar += s.ProgressFilled.Render("█")
		} else {
			bar += s.ProgressEmpty.Render("░")
		}
	}

	return bar
}

// CompressionLabel returns a styled compression label
func (s Styles) CompressionLabel(compression int) string {
	switch compression {
	case 1:
		return s.Value.Render("Real-time")
	case 10:
		return s.Value.Foreground(s.Warning).Render("10x Compressed")
	case 100:
		return s.Value.Foreground(s.Warning).Render("100x Compressed")
	case 1000:
		return s.Value.Foreground(s.Error).Render("1000x Compressed")
	default:
		return s.Value.Foreground(s.Muted).Render("Unknown")
	}
}

// StatusLabel returns a styled status label
func (s Styles) StatusLabel(status string) string {
	switch status {
	case "Ready":
		return s.StatusReady.Render(status)
	case "Running", "Initialization", "Cluster Generation", "Sequencing", "Base Calling":
		return s.StatusRunning.Render(status)
	case "Complete":
		return s.StatusComplete.Render(status)
	case "Failed":
		return s.Value.Foreground(s.Error).Render(status)
	default:
		return s.Value.Foreground(s.Muted).Render(status)
	}
}

// SequencerIcon returns an icon/emoji for sequencer types
func (s Styles) SequencerIcon(model string) string {
	switch model {
	case "NovaSeq 6000":
		return "🚀" // High throughput
	case "NextSeq 2000":
		return "⚡" // Mid throughput
	case "MiSeq":
		return "🔬" // Benchtop
	default:
		return "🧬"
	}
}
