package tui

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	tea "github.com/charmbracelet/bubbletea"
	"github.com/scttfrdmn/phantom-science/internal/timing"
	"github.com/scttfrdmn/phantom-science/pkg/types"
)

// testConfig is a small MiSeq run writing into a temporary directory.
func testConfig(t *testing.T) types.SequencerConfig {
	t.Helper()
	return types.SequencerConfig{
		SequencerModel: "MiSeq",
		FlowCellType:   "v3",
		SampleName:     "TestSample",
		Lane:           1,
		ReadLength:     50,
		PairedEnd:      true,
		OutputDir:      t.TempDir(),
		TotalBatches:   2,
		ReadsPerBatch:  10,
		BaseQuality:    30,
		Seed:           7,
	}
}

// sized returns m after a window-size message, so View renders panels.
func sized(m Model) Model {
	updated, _ := m.Update(tea.WindowSizeMsg{Width: 200, Height: 80})
	return updated.(Model)
}

// flat collapses runs of whitespace, so label padding does not matter.
func flat(s string) string { return strings.Join(strings.Fields(s), " ") }

func press(m Model, k tea.KeyMsg) Model {
	updated, _ := m.Update(k)
	return updated.(Model)
}

var (
	keyTab   = tea.KeyMsg{Type: tea.KeyTab}
	keyEnter = tea.KeyMsg{Type: tea.KeyEnter}
	keyUp    = tea.KeyMsg{Type: tea.KeyUp}
	keyDown  = tea.KeyMsg{Type: tea.KeyDown}
	keyLeft  = tea.KeyMsg{Type: tea.KeyLeft}
	keyRight = tea.KeyMsg{Type: tea.KeyRight}
)

func TestView_LoadingBeforeWindowSize(t *testing.T) {
	if got := New(testConfig(t), nil, "1.2.3").View(); got != "Loading..." {
		t.Errorf("View before WindowSizeMsg = %q, want Loading...", got)
	}
}

// The header showed a hard-coded v0.3.0 regardless of the binary's version.
func TestHeader_ShowsPassedVersion(t *testing.T) {
	view := sized(New(testConfig(t), nil, "9.8.7")).View()
	if !strings.Contains(view, "Phantom Science Sequencer Simulator v9.8.7") {
		t.Errorf("header should show the passed version; got:\n%s", view)
	}
	if strings.Contains(view, "v0.3.0") {
		t.Error("header still shows the hard-coded v0.3.0")
	}
}

// The time-compression highlight compared against (i+1)*10, so only 10x
// could ever be highlighted correctly (and 100x/1000x never were).
func TestOutputPanel_HighlightsSelectedCompression(t *testing.T) {
	for _, selected := range compressionOptions {
		t.Run(fmt.Sprintf("%dx", selected), func(t *testing.T) {
			m := New(testConfig(t), nil, "")
			m.compression = selected
			panel := m.renderOutputPanel()
			for _, c := range compressionOptions {
				marked := strings.Contains(panel, fmt.Sprintf("[%dx]", c))
				if c == selected && !marked {
					t.Errorf("selected %dx not marked; panel:\n%s", c, panel)
				}
				if c != selected && marked {
					t.Errorf("%dx marked although %dx is selected; panel:\n%s", c, selected, panel)
				}
			}
		})
	}
}

func TestOutputPanel_ArrowsChangeCompression(t *testing.T) {
	m := New(testConfig(t), nil, "")
	m.currentSection = SectionOutput
	if m.compression != timing.Compressed100x {
		t.Fatalf("default compression = %d, want 100", m.compression)
	}
	m = press(m, keyRight)
	if m.compression != timing.Compressed1000x {
		t.Errorf("right from 100x = %d, want 1000", m.compression)
	}
	m = press(m, keyRight) // clamps
	if m.compression != timing.Compressed1000x {
		t.Errorf("right at 1000x = %d, want 1000", m.compression)
	}
	for range 4 {
		m = press(m, keyLeft)
	}
	if m.compression != timing.RealTime {
		t.Errorf("left to the start = %d, want 1", m.compression)
	}
	if sim := timing.NewRunSimulator(m.sequencerConfig, m.compression); sim.Compression != timing.RealTime {
		t.Error("selected compression is not what a run would use")
	}
}

// New used to force the MiSeq choice and updateConfigFromUI then overwrote
// --config's sequencer and flow cell, and replaced total_batches with the
// number of demo sample-sheet rows.
func TestUpdateConfigFromUI_KeepsConfigValues(t *testing.T) {
	for _, tc := range []struct{ model, flowCell string }{
		{"NovaSeq 6000", "S4"},
		{"NextSeq 2000", "P2"},
		{"MiSeq", "v2"},
		{"HiSeq 2500", "HighOutput"}, // not a preset
	} {
		t.Run(tc.model, func(t *testing.T) {
			cfg := testConfig(t)
			cfg.SequencerModel, cfg.FlowCellType, cfg.TotalBatches = tc.model, tc.flowCell, 7
			got := New(cfg, nil, "").updateConfigFromUI().sequencerConfig
			if got.SequencerModel != tc.model || got.FlowCellType != tc.flowCell {
				t.Errorf("sequencer = %s/%s, want %s/%s", got.SequencerModel, got.FlowCellType, tc.model, tc.flowCell)
			}
			if got.TotalBatches != 7 {
				t.Errorf("TotalBatches = %d, want 7 from the config", got.TotalBatches)
			}
			if got.SampleName != cfg.SampleName {
				t.Errorf("SampleName = %q, want %q from the config", got.SampleName, cfg.SampleName)
			}
		})
	}
}

func TestStartRun_UsesConfigValues(t *testing.T) {
	cfg := testConfig(t)
	cfg.SequencerModel, cfg.FlowCellType, cfg.TotalBatches = "NovaSeq 6000", "S4", 3
	m := New(cfg, nil, "")
	m.currentSection = SectionActions
	m = press(m, keyEnter)
	defer m.cancel()
	if !m.runInProgress {
		t.Fatalf("run did not start; err = %v", m.err)
	}
	sc := m.simulator.Config
	if sc.SequencerModel != "NovaSeq 6000" || sc.FlowCellType != "S4" || sc.TotalBatches != 3 {
		t.Errorf("simulator config = %s/%s/%d batches, want NovaSeq 6000/S4/3", sc.SequencerModel, sc.FlowCellType, sc.TotalBatches)
	}
	if m.progress.TotalFiles != 6 {
		t.Errorf("TotalFiles = %d, want 6 (3 paired-end batches)", m.progress.TotalFiles)
	}
}

func TestSequencerSelection_AppliesPresetOnlyWhenChosen(t *testing.T) {
	cfg := testConfig(t)
	cfg.FlowCellType = "v2"
	m := New(cfg, nil, "")
	if m.sequencerChoice != ChoiceMiSeq {
		t.Fatalf("choice = %d, want MiSeq from the config", m.sequencerChoice)
	}
	if !strings.Contains(m.renderSequencerPanel(), "MiSeq (v2)") {
		t.Errorf("selected preset should show the configured flow cell; got:\n%s", m.renderSequencerPanel())
	}

	m = press(m, keyUp)
	if m.sequencerChoice != ChoiceNextSeq {
		t.Fatalf("up from MiSeq = %d, want NextSeq", m.sequencerChoice)
	}
	got := m.updateConfigFromUI().sequencerConfig
	if got.SequencerModel != "NextSeq 2000" || got.FlowCellType != "P3" {
		t.Errorf("after choosing NextSeq: %s/%s, want NextSeq 2000/P3", got.SequencerModel, got.FlowCellType)
	}
	m = press(press(m, keyUp), keyUp) // clamps at NovaSeq
	if m.sequencerChoice != ChoiceNovaSeq {
		t.Errorf("choice = %d, want NovaSeq", m.sequencerChoice)
	}
	m = press(press(press(m, keyDown), keyDown), keyDown) // clamps at MiSeq
	if m.sequencerChoice != ChoiceMiSeq {
		t.Errorf("choice = %d, want MiSeq", m.sequencerChoice)
	}
	if got := m.updateConfigFromUI().sequencerConfig.FlowCellType; got != "v3" {
		t.Errorf("re-chosen MiSeq flow cell = %s, want the preset v3", got)
	}
}

func TestSequencerPanel_ShowsNonPresetConfig(t *testing.T) {
	cfg := testConfig(t)
	cfg.SequencerModel, cfg.FlowCellType = "HiSeq 2500", "Rapid"
	m := New(cfg, nil, "")
	if m.sequencerChoice != ChoiceConfig {
		t.Fatalf("choice = %d, want ChoiceConfig", m.sequencerChoice)
	}
	if panel := m.renderSequencerPanel(); !strings.Contains(panel, "● HiSeq 2500 (Rapid)") {
		t.Errorf("panel should show the configured sequencer as selected; got:\n%s", panel)
	}
}

func TestSequencerAndCompression_LockedDuringRun(t *testing.T) {
	m := New(testConfig(t), nil, "")
	m.runInProgress = true
	m = press(m, keyUp)
	m.currentSection = SectionOutput
	m = press(m, keyLeft)
	if m.sequencerChoice != ChoiceMiSeq || m.sequencerChanged || m.compression != timing.Compressed100x {
		t.Error("arrow keys changed settings during a run, which the running simulation would not pick up")
	}
}

// Enter in the Samples panel used to enter an editing mode that swallowed
// Tab until Esc. Tab must reach every section from every section.
func TestTab_CyclesAllSectionsFromEverySection(t *testing.T) {
	for start := range numSections {
		m := New(testConfig(t), nil, "")
		m.currentSection = start
		m = press(m, keyEnter) // must not trap focus (Actions would start a run)
		if m.cancel != nil {
			m.cancel()
		}
		m.runInProgress = false
		seen := map[int]bool{}
		for range numSections {
			m = press(m, keyTab)
			seen[m.currentSection] = true
		}
		if len(seen) != numSections || m.currentSection != start {
			t.Errorf("from section %d: Tab visited %v and ended at %d", start, seen, m.currentSection)
		}
		for range numSections {
			m = press(m, tea.KeyMsg{Type: tea.KeyShiftTab})
		}
		if m.currentSection != start {
			t.Errorf("from section %d: Shift+Tab cycle ended at %d", start, m.currentSection)
		}
	}
}

func TestQuitKeys(t *testing.T) {
	m := New(testConfig(t), nil, "")
	for _, k := range []tea.KeyMsg{{Type: tea.KeyCtrlC}, {Type: tea.KeyRunes, Runes: []rune{'q'}}} {
		_, cmd := m.Update(k)
		if cmd == nil {
			t.Fatalf("%s returned no command", k)
		}
		if _, ok := cmd().(tea.QuitMsg); !ok {
			t.Errorf("%s did not quit", k)
		}
	}
}

func TestDisplay_ReflectsCompressionAndPairedEnd(t *testing.T) {
	cfg := testConfig(t)
	cfg.Output.Compression = "gzip"
	cfg.PairedEnd = false
	view := flat(sized(New(cfg, nil, "")).View())
	for _, want := range []string{"Compression: [gzip]", "Paired-End: [✗]", "TestSample_S1_L001_R1_001.fastq.gz", "2 FASTQ files"} {
		if !strings.Contains(view, want) {
			t.Errorf("view missing %q; got:\n%s", want, view)
		}
	}

	cfg.Output.Compression = ""
	cfg.PairedEnd = true
	view = flat(sized(New(cfg, nil, "")).View())
	for _, want := range []string{"Compression: [none]", "Paired-End: [✓]", "TestSample_S1_L001_R1_001.fastq", "4 FASTQ files"} {
		if !strings.Contains(view, want) {
			t.Errorf("view missing %q; got:\n%s", want, view)
		}
	}
	for _, gone := range []string{"gzip", ".fastq.gz", "v1.5", "PhiX", "Real-time", "Load Config", "Save Config", "[Help]", "Add Sample", "Import CSV", "Strategy"} {
		if strings.Contains(view, gone) {
			t.Errorf("view still shows %q", gone)
		}
	}
}

func TestPreview_ReflectsChosenSequencer(t *testing.T) {
	m := New(testConfig(t), nil, "")
	m = press(m, keyUp) // NextSeq 2000 P3
	if p := m.renderPreview(); !strings.Contains(p, "NextSeq 2000 (P3)") {
		t.Errorf("preview should describe the chosen sequencer; got:\n%s", p)
	}
}

// The panel showed only the phantom config's enable flag, so fingerprinting
// set in the sequencer config showed as off although the run applied it.
func TestFingerprintingPanel_ShowsResolvedSettings(t *testing.T) {
	cfg := testConfig(t)
	cfg.Fingerprinting = types.FingerprintConfig{Enable: true, SteganographyLayers: []string{"quality_scores"}}
	cfg.Enterprise.CustomerID = "ACME"
	panel := flat(New(cfg, nil, "").renderFingerprintingPanel())
	for _, want := range []string{"Enabled: [✓]", "[ACME]", "[quality_scores]"} {
		if !strings.Contains(panel, want) {
			t.Errorf("panel missing %q; got:\n%s", want, panel)
		}
	}

	phantom := &types.PhantomConfig{Enterprise: types.EnterpriseConfig{
		CustomerID:     "PHANTOM_CUST",
		Fingerprinting: types.FingerprintConfig{Enable: true},
	}}
	m := New(testConfig(t), phantom, "")
	if !m.fingerprintingEnabled || m.customerID != "PHANTOM_CUST" {
		t.Errorf("phantom config: enabled=%v customer=%q", m.fingerprintingEnabled, m.customerID)
	}

	if m := New(testConfig(t), nil, ""); m.fingerprintingEnabled {
		t.Error("fingerprinting shown as enabled with no config enabling it")
	}
}

func TestSampleSheet_ShowsConfiguredSample(t *testing.T) {
	cfg := testConfig(t)
	cfg.SampleName, cfg.Lane = "Real_Sample", 3
	sheet := flat(New(cfg, nil, "").renderSampleSheet())
	if !strings.Contains(sheet, "Real_Sample 3 — —") {
		t.Errorf("sample sheet should show sample_name and lane; got:\n%s", sheet)
	}
	if strings.Contains(sheet, "Patient_001") {
		t.Error("sample sheet still shows demo samples")
	}
}

// A file that fails to generate was printed with fmt.Printf into the
// alternate screen; it must appear in View.
func TestRun_GenerationErrorShownInView(t *testing.T) {
	cfg := testConfig(t)
	cfg.SampleName = "../escape" // file path leaves output_dir and is rejected
	m := sized(New(cfg, nil, ""))
	m.compression = 1_000_000
	m.currentSection = SectionActions
	m, cmd := m.handleActionSelection()
	if !m.runInProgress || cmd == nil {
		t.Fatalf("run did not start; err = %v", m.err)
	}
	updated, _ := m.Update(cmd()) // runs the simulation to completion
	m = updated.(Model)

	view := m.View()
	if !strings.Contains(view, "Error:") || !strings.Contains(view, "file(s) failed") {
		t.Errorf("view should show the generation error; got:\n%s", view)
	}
	if m.progress.FilesGenerated != 0 {
		t.Errorf("FilesGenerated = %d, want 0", m.progress.FilesGenerated)
	}
}

func TestRun_OutputDirErrorShownInView(t *testing.T) {
	cfg := testConfig(t)
	file := filepath.Join(cfg.OutputDir, "not-a-dir")
	if err := os.WriteFile(file, nil, 0600); err != nil {
		t.Fatal(err)
	}
	cfg.OutputDir = filepath.Join(file, "out")
	m := sized(New(cfg, nil, ""))
	m.currentSection = SectionActions
	m = press(m, keyEnter)
	if m.runInProgress {
		t.Fatal("run started although the output directory could not be created")
	}
	if view := m.View(); !strings.Contains(view, "creating output directory") {
		t.Errorf("view should show the output directory error; got:\n%s", view)
	}
}

// A run writes total_batches batches and the panels show what was written.
func TestRun_WritesFilesAndReportsProgress(t *testing.T) {
	cfg := testConfig(t)
	cfg.Output.Compression = "gzip"
	m := sized(New(cfg, nil, "0.9.0"))
	m.compression = 1_000_000
	m.currentSection = SectionActions
	m, cmd := m.handleActionSelection()
	if cmd == nil {
		t.Fatalf("run did not start; err = %v", m.err)
	}
	msg := cmd()
	if _, ok := msg.(simulationCompleteMsg); !ok {
		t.Fatalf("simulation ended with %#v", msg)
	}
	updated, _ := m.Update(msg)
	m = updated.(Model)

	if m.err != nil {
		t.Fatalf("run error: %v", m.err)
	}
	files, _ := filepath.Glob(filepath.Join(cfg.OutputDir, "*.fastq.gz"))
	if len(files) != 4 {
		t.Errorf("wrote %d .fastq.gz files (%v), want 4", len(files), files)
	}
	var total int64
	for _, f := range files {
		fi, err := os.Stat(f)
		if err != nil {
			t.Fatal(err)
		}
		total += fi.Size()
	}
	if m.progress.FilesGenerated != 4 || m.progress.TotalFiles != 4 || m.progress.BytesGenerated != total {
		t.Errorf("progress = %d/%d files, %d bytes; want 4/4, %d bytes",
			m.progress.FilesGenerated, m.progress.TotalFiles, m.progress.BytesGenerated, total)
	}
	view := flat(m.View())
	for _, want := range []string{"Status: Complete", "Files: 4/4", ".fastq.gz"} {
		if !strings.Contains(view, want) {
			t.Errorf("view missing %q; got:\n%s", want, view)
		}
	}
}

// Stopping a run and starting another must not let the first run's
// cancellation message end the second.
func TestRun_StaleMessageIgnored(t *testing.T) {
	m := New(testConfig(t), nil, "")
	m.currentSection = SectionActions
	m = press(m, keyEnter) // start run 1
	first := m.runID
	m = press(m, keyEnter) // stop
	if m.runInProgress || m.progress.Stage != stageStopped {
		t.Fatalf("after stop: running=%v stage=%s", m.runInProgress, m.progress.Stage)
	}
	m = press(m, keyEnter) // start run 2
	defer m.cancel()
	updated, _ := m.Update(simulationErrorMsg{runID: first, err: context.Canceled})
	m = updated.(Model)
	if !m.runInProgress {
		t.Error("the stopped run's message ended the new run")
	}
	if m.err != nil {
		t.Errorf("cancellation reported as an error: %v", m.err)
	}
}

func TestRun_SimulatorErrorShownInView(t *testing.T) {
	m := sized(New(testConfig(t), nil, ""))
	m.currentSection = SectionActions
	m = press(m, keyEnter)
	defer m.cancel()
	updated, _ := m.Update(simulationErrorMsg{runID: m.runID, err: fmt.Errorf("simulator broke")})
	m = updated.(Model)
	if m.runInProgress {
		t.Error("run still in progress after an error")
	}
	if view := flat(m.View()); !strings.Contains(view, "Error: simulator broke") || !strings.Contains(view, "Status: Failed") {
		t.Errorf("view should show the failure; got:\n%s", view)
	}
}

func TestParametersPanel_ShowsConfigValues(t *testing.T) {
	cfg := testConfig(t)
	cfg.HMMModelPath = "/very/long/deeply/nested/path/to/models/GBR.hmm.json"
	cfg.Population = "GBR"
	cfg.ReadsPerBatch, cfg.TotalBatches = 123, 9
	panel := New(cfg, nil, "").renderParametersPanel()
	for _, want := range []string{"[GBR.hmm.json]", "[GBR]", "[50 bp]", "[123]", "[9]", "[7]", "[illumina]"} {
		if !strings.Contains(panel, want) {
			t.Errorf("parameters panel missing %q; got:\n%s", want, panel)
		}
	}

	cfg.HMMModelPath, cfg.Population, cfg.Seed = "", "", 0
	panel = New(cfg, nil, "").renderParametersPanel()
	if !strings.Contains(panel, "[random]") || strings.Count(panel, "[—]") != 2 {
		t.Errorf("unset seed/HMM/population should show random and dashes; got:\n%s", panel)
	}
}
