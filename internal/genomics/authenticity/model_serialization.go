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

// adaptiveHMMModelData is the serializable snapshot of AdaptiveHMMModel.
// Excludes: mu, rng, complexityMemory (runtime sliding window),
// currentComplexity, currentGCContent (runtime metrics).
type adaptiveHMMModelData struct {
	K                    int                `json:"k"`
	StateSpace           []string           `json:"state_space"`
	TransitionMatrix     [][]float64        `json:"transition_matrix"`
	EmissionMatrix       [][]float64        `json:"emission_matrix"`
	StartProbabilities   []float64          `json:"start_probabilities"`
	TotalKmers           uint64             `json:"total_kmers"`
	Trained              bool               `json:"trained"`
	ComplexityWindowSize int                `json:"complexity_window_size"`
	TargetComplexity     float64            `json:"target_complexity"`
	GCContentTarget      float64            `json:"gc_content_target"`
	GCContentTolerance   float64            `json:"gc_content_tolerance"`
	HomopolymerMaxLength int                `json:"homopolymer_max_length"`
	CpGDepletionRatio    float64            `json:"cpg_depletion_ratio"`
	AdaptiveWeights      map[string]float64 `json:"adaptive_weights"`
}

const (
	envelopeFormat  = "phantom-science-hmm"
	envelopeVersion = "1"
)

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
func LoadHMMKmerModel(path string) (*HMMKmerModel, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("hmm_kmer_model: read file: %w", err)
	}

	var env ModelEnvelope
	if err := json.Unmarshal(raw, &env); err != nil {
		return nil, fmt.Errorf("hmm_kmer_model: unmarshal envelope: %w", err)
	}

	if env.Format != envelopeFormat {
		return nil, fmt.Errorf("hmm_kmer_model: unexpected format %q", env.Format)
	}
	if env.Version != envelopeVersion {
		return nil, fmt.Errorf("hmm_kmer_model: unsupported version %q", env.Version)
	}
	if env.ModelType != "HMMKmerModel" {
		return nil, fmt.Errorf("hmm_kmer_model: wrong model type %q", env.ModelType)
	}

	var data hmmKmerModelData
	if err := json.Unmarshal(env.Data, &data); err != nil {
		return nil, fmt.Errorf("hmm_kmer_model: unmarshal data: %w", err)
	}

	m := NewHMMKmerModel(data.K)
	if m == nil {
		return nil, fmt.Errorf("hmm_kmer_model: invalid k=%d in saved model", data.K)
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

// Save serializes the trained AdaptiveHMMModel to a JSON file at path.
// Returns an error if the model has not been trained.
// The recommended file extension is .hmm.json.
func (m *AdaptiveHMMModel) Save(path string) error {
	m.mu.RLock()
	defer m.mu.RUnlock()

	if !m.trained {
		return fmt.Errorf("adaptive_hmm_model: cannot save untrained model")
	}

	data := adaptiveHMMModelData{
		K:                    m.k,
		StateSpace:           m.stateSpace,
		TransitionMatrix:     m.transitionMatrix,
		EmissionMatrix:       m.emissionMatrix,
		StartProbabilities:   m.startProbabilities,
		TotalKmers:           m.totalKmers,
		Trained:              m.trained,
		ComplexityWindowSize: m.complexityWindowSize,
		TargetComplexity:     m.targetComplexity,
		GCContentTarget:      m.gcContentTarget,
		GCContentTolerance:   m.gcContentTolerance,
		HomopolymerMaxLength: m.homopolymerMaxLength,
		CpGDepletionRatio:    m.cpgDepletionRatio,
		AdaptiveWeights:      m.adaptiveWeights,
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

// LoadAdaptiveHMMModel deserializes an AdaptiveHMMModel from the JSON file at path.
func LoadAdaptiveHMMModel(path string) (*AdaptiveHMMModel, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("adaptive_hmm_model: read file: %w", err)
	}

	var env ModelEnvelope
	if err := json.Unmarshal(raw, &env); err != nil {
		return nil, fmt.Errorf("adaptive_hmm_model: unmarshal envelope: %w", err)
	}

	if env.Format != envelopeFormat {
		return nil, fmt.Errorf("adaptive_hmm_model: unexpected format %q", env.Format)
	}
	if env.Version != envelopeVersion {
		return nil, fmt.Errorf("adaptive_hmm_model: unsupported version %q", env.Version)
	}
	if env.ModelType != "AdaptiveHMMModel" {
		return nil, fmt.Errorf("adaptive_hmm_model: wrong model type %q", env.ModelType)
	}

	var data adaptiveHMMModelData
	if err := json.Unmarshal(env.Data, &data); err != nil {
		return nil, fmt.Errorf("adaptive_hmm_model: unmarshal data: %w", err)
	}

	m := NewAdaptiveHMMModel(data.K)
	if m == nil {
		return nil, fmt.Errorf("adaptive_hmm_model: invalid k=%d in saved model", data.K)
	}

	m.stateSpace = data.StateSpace
	m.transitionMatrix = data.TransitionMatrix
	m.emissionMatrix = data.EmissionMatrix
	m.startProbabilities = data.StartProbabilities
	m.totalKmers = data.TotalKmers
	m.trained = data.Trained
	m.complexityWindowSize = data.ComplexityWindowSize
	m.targetComplexity = data.TargetComplexity
	m.gcContentTarget = data.GCContentTarget
	m.gcContentTolerance = data.GCContentTolerance
	m.homopolymerMaxLength = data.HomopolymerMaxLength
	m.cpgDepletionRatio = data.CpGDepletionRatio
	m.adaptiveWeights = data.AdaptiveWeights
	m.complexityMemory = make([]string, 0)
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
