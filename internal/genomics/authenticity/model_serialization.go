package authenticity

import (
	"encoding/json"
	"fmt"
	"math/rand/v2"
	"os"
	"time"
)

// SequenceGenerator is the minimal interface satisfied by both HMMKmerModel
// and AdaptiveHMMModel, used by the genomics engine.
type SequenceGenerator interface {
	GenerateSequence(length int) (*GenerationResult, error)
	IsTrained() bool
}

// ModelEnvelope is the top-level JSON container for any serialized model.
type ModelEnvelope struct {
	Format       string          `json:"format"`
	Version      string          `json:"version"`
	ModelType    string          `json:"model_type"`
	SerializedAt time.Time       `json:"serialized_at"`
	Population   string          `json:"population,omitempty"`
	Data         json.RawMessage `json:"data"`
}

// hmmKmerModelData is the serializable snapshot of HMMKmerModel.
// Excludes: mu (sync.RWMutex), rng (*rand.Rand) — both recreated on load.
type hmmKmerModelData struct {
	K                  int         `json:"k"`
	StateSpace         []string    `json:"state_space"`
	TransitionMatrix   [][]float64 `json:"transition_matrix"`
	EmissionMatrix     [][]float64 `json:"emission_matrix"`
	StartProbabilities []float64   `json:"start_probabilities"`
	TotalKmers         uint64      `json:"total_kmers"`
	Trained            bool        `json:"trained"`
}

// adaptiveHMMModelData is the serializable snapshot of AdaptiveHMMModel
// (format version 2). It holds the raw training counts, so a loaded model can
// be trained further, and the finalized probabilities estimated from them
// with PseudoCount. Excludes: mu, rng, complexityMemory and generation
// statistics (runtime state).
type adaptiveHMMModelData struct {
	K                    int         `json:"k"`
	StateSpace           []string    `json:"state_space"`
	PseudoCount          float64     `json:"pseudocount"`
	TransitionCounts     [][]float64 `json:"transition_counts"`
	TransitionMatrix     [][]float64 `json:"transition_matrix"`
	StartProbabilities   []float64   `json:"start_probabilities"`
	TotalKmers           uint64      `json:"total_kmers"`
	Trained              bool        `json:"trained"`
	ComplexityWindowSize int         `json:"complexity_window_size"`
	TargetComplexity     float64     `json:"target_complexity"`
}

// Envelope format and versions. Version 2 (this release) changed how models
// are trained: transition and start probabilities use pseudocount
// PseudoCount, and adaptive models store their raw counts. Version 1
// HMMKmerModel files are still loaded (their probabilities are normalised);
// version 1 AdaptiveHMMModel files are rejected because they were trained
// with the wrong context and no normalisation (generation from them was
// degenerate, e.g. all C), and must be retrained.
const (
	envelopeFormat  = "phantom-science-hmm"
	envelopeVersion = "2"
)

// readEnvelope reads and checks the envelope of the model file at path.
func readEnvelope(component, path, modelType string) (*ModelEnvelope, error) {
	raw, err := os.ReadFile(path) // #nosec G304 -- model path chosen by the user
	if err != nil {
		return nil, fmt.Errorf("%s: read file: %w", component, err)
	}
	var env ModelEnvelope
	if err := json.Unmarshal(raw, &env); err != nil {
		return nil, fmt.Errorf("%s: unmarshal envelope: %w", component, err)
	}
	if env.Format != envelopeFormat {
		return nil, fmt.Errorf("%s: unexpected format %q", component, env.Format)
	}
	if env.ModelType != modelType {
		return nil, fmt.Errorf("%s: wrong model type %q", component, env.ModelType)
	}
	return &env, nil
}

// Save serializes the trained HMMKmerModel to a JSON file at path.
// Returns an error if the model has not been trained.
// The recommended file extension is .hmm.json.
func (m *HMMKmerModel) Save(path string) error {
	m.mu.RLock()
	defer m.mu.RUnlock()

	if !m.trained {
		return fmt.Errorf("hmm_kmer_model: cannot save untrained model")
	}

	data := hmmKmerModelData{
		K:                  m.k,
		StateSpace:         m.stateSpace,
		TransitionMatrix:   m.transitionMatrix,
		EmissionMatrix:     m.emissionMatrix,
		StartProbabilities: m.startProbabilities,
		TotalKmers:         m.totalKmers,
		Trained:            m.trained,
	}

	raw, err := json.Marshal(data)
	if err != nil {
		return fmt.Errorf("hmm_kmer_model: marshal data: %w", err)
	}

	env := ModelEnvelope{
		Format:       envelopeFormat,
		Version:      envelopeVersion,
		ModelType:    "HMMKmerModel",
		SerializedAt: time.Now().UTC(),
		Data:         json.RawMessage(raw),
	}

	out, err := json.MarshalIndent(env, "", "  ")
	if err != nil {
		return fmt.Errorf("hmm_kmer_model: marshal envelope: %w", err)
	}

	if err := os.WriteFile(path, out, 0600); err != nil {
		return fmt.Errorf("hmm_kmer_model: write file: %w", err)
	}

	return nil
}

// LoadHMMKmerModel deserializes a HMMKmerModel from the JSON file at path.
// Format versions 1 and 2 are accepted.
func LoadHMMKmerModel(path string) (*HMMKmerModel, error) {
	env, err := readEnvelope("hmm_kmer_model", path, "HMMKmerModel")
	if err != nil {
		return nil, err
	}
	if env.Version != "1" && env.Version != envelopeVersion {
		return nil, fmt.Errorf("hmm_kmer_model: unsupported version %q", env.Version)
	}

	var data hmmKmerModelData
	if err := json.Unmarshal(env.Data, &data); err != nil {
		return nil, fmt.Errorf("hmm_kmer_model: unmarshal data: %w", err)
	}
	if data.K < 2 || data.K > 8 {
		return nil, fmt.Errorf("hmm_kmer_model: invalid k=%d in saved model", data.K)
	}

	m := NewHMMKmerModel(data.K)
	if m == nil {
		return nil, fmt.Errorf("hmm_kmer_model: invalid k=%d in saved model", data.K)
	}
	if err := validateChain("hmm_kmer_model", data.K, data.StateSpace, nil, data.TransitionMatrix, data.StartProbabilities, false); err != nil {
		return nil, fmt.Errorf("hmm_kmer_model: malformed model: %w", err)
	}
	if len(data.EmissionMatrix) != len(data.StateSpace) {
		return nil, fmt.Errorf("hmm_kmer_model: malformed model: %d emission rows for %d states", len(data.EmissionMatrix), len(data.StateSpace))
	}
	for i, row := range data.EmissionMatrix {
		if len(row) != 42 {
			return nil, fmt.Errorf("hmm_kmer_model: malformed model: emission row %d has %d values, want 42", i, len(row))
		}
	}

	m.stateSpace = data.StateSpace
	m.transitionMatrix = data.TransitionMatrix
	m.emissionMatrix = data.EmissionMatrix
	m.startProbabilities = data.StartProbabilities
	m.totalKmers = data.TotalKmers
	m.trained = data.Trained
	m.rng = rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)) // #nosec G404

	if err := m.ValidateModel(); err != nil {
		return nil, fmt.Errorf("hmm_kmer_model: validation after load: %w", err)
	}

	return m, nil
}

// Save serializes the trained AdaptiveHMMModel to a JSON file at path,
// finalizing it first. Returns an error if the model has not been trained.
// The recommended file extension is .hmm.json.
func (m *AdaptiveHMMModel) Save(path string) error {
	m.mu.Lock()
	defer m.mu.Unlock()

	if !m.trained {
		return fmt.Errorf("adaptive_hmm_model: cannot save untrained model")
	}
	m.finalizeLocked()

	data := adaptiveHMMModelData{
		K:                    m.k,
		StateSpace:           m.stateSpace,
		PseudoCount:          PseudoCount,
		TransitionCounts:     m.transitionCounts,
		TransitionMatrix:     m.transitionMatrix,
		StartProbabilities:   m.startProbabilities,
		TotalKmers:           m.totalKmers,
		Trained:              m.trained,
		ComplexityWindowSize: m.complexityWindowSize,
		TargetComplexity:     m.targetComplexity,
	}

	raw, err := json.Marshal(data)
	if err != nil {
		return fmt.Errorf("adaptive_hmm_model: marshal data: %w", err)
	}

	env := ModelEnvelope{
		Format:       envelopeFormat,
		Version:      envelopeVersion,
		ModelType:    "AdaptiveHMMModel",
		SerializedAt: time.Now().UTC(),
		Data:         json.RawMessage(raw),
	}

	out, err := json.MarshalIndent(env, "", "  ")
	if err != nil {
		return fmt.Errorf("adaptive_hmm_model: marshal envelope: %w", err)
	}

	if err := os.WriteFile(path, out, 0600); err != nil {
		return fmt.Errorf("adaptive_hmm_model: write file: %w", err)
	}

	return nil
}

// LoadAdaptiveHMMModel deserializes an AdaptiveHMMModel from the JSON file at
// path. Only format version 2 is accepted; see envelopeVersion.
func LoadAdaptiveHMMModel(path string) (*AdaptiveHMMModel, error) {
	env, err := readEnvelope("adaptive_hmm_model", path, "AdaptiveHMMModel")
	if err != nil {
		return nil, err
	}
	switch env.Version {
	case envelopeVersion:
	case "1":
		return nil, fmt.Errorf("adaptive_hmm_model: %s is a version 1 model; version 1 adaptive models were "+
			"trained with an off-by-one context and no smoothing (their output was degenerate), "+
			"so they are no longer loaded: retrain it with `phantom-science train`", path)
	default:
		return nil, fmt.Errorf("adaptive_hmm_model: unsupported version %q", env.Version)
	}

	var data adaptiveHMMModelData
	if err := json.Unmarshal(env.Data, &data); err != nil {
		return nil, fmt.Errorf("adaptive_hmm_model: unmarshal data: %w", err)
	}
	if data.K < 2 || data.K > 8 {
		return nil, fmt.Errorf("adaptive_hmm_model: invalid k=%d in saved model", data.K)
	}
	if data.PseudoCount != PseudoCount {
		return nil, fmt.Errorf("adaptive_hmm_model: model pseudocount %v, this version uses %v", data.PseudoCount, PseudoCount)
	}
	if err := validateChain("adaptive_hmm_model", data.K, data.StateSpace, data.TransitionCounts, data.TransitionMatrix, data.StartProbabilities, true); err != nil {
		return nil, fmt.Errorf("adaptive_hmm_model: validation after load: %w", err)
	}

	m := NewAdaptiveHMMModel(data.K)
	if m == nil {
		return nil, fmt.Errorf("adaptive_hmm_model: invalid k=%d in saved model", data.K)
	}
	m.stateSpace = data.StateSpace
	m.transitionCounts = data.TransitionCounts
	m.transitionMatrix = data.TransitionMatrix
	m.startProbabilities = data.StartProbabilities
	m.finalized = true
	m.totalKmers = data.TotalKmers
	m.trained = data.Trained
	if data.ComplexityWindowSize > 0 {
		m.complexityWindowSize = data.ComplexityWindowSize
	}
	m.targetComplexity = data.TargetComplexity
	m.rng = rand.New(rand.NewPCG(uint64(time.Now().UnixNano()), 0)) // #nosec G404

	if err := m.ValidateModel(); err != nil {
		return nil, fmt.Errorf("adaptive_hmm_model: validation after load: %w", err)
	}

	return m, nil
}

// LoadModel loads any supported HMM model from path and returns a SequenceGenerator.
// Peeks at ModelEnvelope.ModelType to dispatch to the correct loader.
func LoadModel(path string) (SequenceGenerator, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("load_model: read %s: %w", path, err)
	}
	var env ModelEnvelope
	if err := json.Unmarshal(raw, &env); err != nil {
		return nil, fmt.Errorf("load_model: parse envelope: %w", err)
	}
	switch env.ModelType {
	case "AdaptiveHMMModel":
		return LoadAdaptiveHMMModel(path)
	case "HMMKmerModel":
		return LoadHMMKmerModel(path)
	default:
		return nil, fmt.Errorf("load_model: unknown model type %q", env.ModelType)
	}
}
