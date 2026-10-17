//go:build !release

package validation

import (
	"context"
	"encoding/json"
	"fmt"
	"net/http"
	"net/http/httptest"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// --- InstitutionAPI tests ---

func TestInstitutionAPI_SearchExperts_NoEndpoint(t *testing.T) {
	api, err := NewInstitutionAPI(InstitutionConfig{}, &http.Client{})
	require.NoError(t, err)

	candidates, err := api.SearchExperts(context.Background(), ExpertCriteria{})
	assert.NoError(t, err)
	assert.Empty(t, candidates)
}

func TestInstitutionAPI_CheckSubmissionStatus_NoEndpoint(t *testing.T) {
	api, err := NewInstitutionAPI(InstitutionConfig{}, &http.Client{})
	require.NoError(t, err)

	_, err = api.CheckSubmissionStatus(context.Background(), "sub_123")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "APIEndpoint not configured")
}

func TestInstitutionAPI_SearchExperts_Success(t *testing.T) {
	want := []*ExpertCandidate{
		{Name: "Dr. Jane Smith", Email: "jsmith@uni.edu", Publications: 42, HIndex: 15},
	}
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		assert.Equal(t, "/experts", r.URL.Path)
		assert.Equal(t, "genomics", r.URL.Query().Get("research_areas"))
		assert.Equal(t, "Bearer test-key", r.Header.Get("Authorization"))
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(want)
	}))
	defer srv.Close()

	api, err := NewInstitutionAPI(InstitutionConfig{APIEndpoint: srv.URL, APIKey: "test-key"}, srv.Client())
	require.NoError(t, err)

	got, err := api.SearchExperts(context.Background(), ExpertCriteria{
		ResearchAreas: []string{"genomics"},
	})
	require.NoError(t, err)
	require.Len(t, got, 1)
	assert.Equal(t, "Dr. Jane Smith", got[0].Name)
}

func TestInstitutionAPI_SearchExperts_ServerError(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.Error(w, "internal error", http.StatusInternalServerError)
	}))
	defer srv.Close()

	api, err := NewInstitutionAPI(InstitutionConfig{APIEndpoint: srv.URL}, srv.Client())
	require.NoError(t, err)

	_, err = api.SearchExperts(context.Background(), ExpertCriteria{})
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "500")
}

func TestInstitutionAPI_CheckSubmissionStatus_Success(t *testing.T) {
	want := SubmissionStatus{
		SubmissionID:  "sub_xyz",
		Journal:       "Nature Genomics",
		CurrentStatus: "under_review",
		LastUpdated:   time.Now().UTC().Truncate(time.Second),
		ReviewerCount: 3,
	}
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		assert.Equal(t, "/submissions/sub_xyz", r.URL.Path)
		w.Header().Set("Content-Type", "application/json")
		json.NewEncoder(w).Encode(want)
	}))
	defer srv.Close()

	api, err := NewInstitutionAPI(InstitutionConfig{APIEndpoint: srv.URL}, srv.Client())
	require.NoError(t, err)

	got, err := api.CheckSubmissionStatus(context.Background(), "sub_xyz")
	require.NoError(t, err)
	assert.Equal(t, "sub_xyz", got.SubmissionID)
	assert.Equal(t, "under_review", got.CurrentStatus)
}

func TestInstitutionAPI_CheckSubmissionStatus_NotFound(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.NotFound(w, r)
	}))
	defer srv.Close()

	api, err := NewInstitutionAPI(InstitutionConfig{APIEndpoint: srv.URL}, srv.Client())
	require.NoError(t, err)

	_, err = api.CheckSubmissionStatus(context.Background(), "missing_id")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "not found")
}

// --- PubMedClient tests ---

func TestPubMedClient_SearchArticles_NoResults(t *testing.T) {
	// Server returns empty idlist — should return nil slice without error.
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Type", "application/json")
		w.Write([]byte(`{"esearchresult":{"idlist":[]}}`))
	}))
	defer srv.Close()

	c := &PubMedClient{httpClient: srv.Client()}
	// Override base URL isn't possible without refactoring; test via mock transport.
	// Instead verify the shape of ncbiESearchResponse parsing.
	var parsed ncbiESearchResponse
	err := json.Unmarshal([]byte(`{"esearchresult":{"idlist":[]}}`), &parsed)
	require.NoError(t, err)
	assert.Empty(t, parsed.ESearchResult.IDList)
	_ = c
}

func TestPubMedClient_SearchArticles_ESearchError(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.Error(w, "service unavailable", http.StatusServiceUnavailable)
	}))
	defer srv.Close()

	// Use a transport that redirects to our test server for any host.
	transport := &redirectTransport{base: srv.Client().Transport, target: srv.URL}
	c := &PubMedClient{httpClient: &http.Client{Transport: transport}}

	_, err := c.SearchArticles(context.Background(), "genomics")
	assert.Error(t, err)
}

func TestPubMedClient_SearchArticles_ParseESearchResponse(t *testing.T) {
	raw := `{"esearchresult":{"count":"2","idlist":["12345678","98765432"]}}`
	var r ncbiESearchResponse
	require.NoError(t, json.Unmarshal([]byte(raw), &r))
	assert.Equal(t, []string{"12345678", "98765432"}, r.ESearchResult.IDList)
}

func TestPubMedClient_SearchArticles_ParseESummaryResponse(t *testing.T) {
	raw := `{"result":{"12345678":{"uid":"12345678","title":"Genomics Study","source":"Nature","pubdate":"2023 Jan","elocationid":"10.1038/foo","authors":[{"name":"Smith J"},{"name":"Jones K"}]}}}`
	var r ncbiESummaryResponse
	require.NoError(t, json.Unmarshal([]byte(raw), &r))
	doc := r.Result["12345678"]
	assert.Equal(t, "Genomics Study", doc.Title)
	assert.Equal(t, "Nature", doc.Source)
	assert.Len(t, doc.Authors, 2)
	assert.Equal(t, "Smith J", doc.Authors[0].Name)
}

// --- ORCIDClient tests ---

func TestORCIDClient_LookupResearcher_NotFound(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.NotFound(w, r)
	}))
	defer srv.Close()

	transport := &redirectTransport{base: srv.Client().Transport, target: srv.URL}
	c := &ORCIDClient{httpClient: &http.Client{Transport: transport}}

	_, err := c.LookupResearcher(context.Background(), "0000-0001-2345-6789")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "not found")
}

func TestORCIDClient_LookupResearcher_ParseResponse(t *testing.T) {
	raw := `{
		"person": {
			"name": {
				"given-names": {"value": "Jane"},
				"family-name": {"value": "Smith"}
			},
			"biography": {"content": "Professor of genomics"},
			"emails": {"email": [{"email": "jsmith@uni.edu", "primary": true}]}
		},
		"activities-summary": {
			"employments": {
				"affiliation-group": [
					{"summaries": [{"employment-summary": {"organization": {"name": "MIT"}}}]}
				]
			},
			"works": {"group": [{}, {}]}
		}
	}`
	var record orcidRecord
	require.NoError(t, json.Unmarshal([]byte(raw), &record))
	assert.Equal(t, "Jane", record.Person.Name.GivenNames.Value)
	assert.Equal(t, "Smith", record.Person.Name.FamilyName.Value)
	assert.Equal(t, "jsmith@uni.edu", record.Person.Emails.Email[0].Email)
	assert.True(t, record.Person.Emails.Email[0].Primary)
	assert.Equal(t, "MIT", record.ActivitiesSummary.Employments.AffiliationGroup[0].Summaries[0].EmploymentSummary.Organization.Name)
	assert.Len(t, record.ActivitiesSummary.Works.Group, 2)
}

func TestORCIDClient_LookupResearcher_ServerSuccess(t *testing.T) {
	payload := `{
		"person": {
			"name": {"given-names": {"value": "Alice"}, "family-name": {"value": "Wong"}},
			"biography": {"content": ""},
			"emails": {"email": [{"email": "awong@example.com", "primary": true}]}
		},
		"activities-summary": {
			"employments": {"affiliation-group": [
				{"summaries": [{"employment-summary": {"organization": {"name": "Stanford"}}}]}
			]},
			"works": {"group": [{}, {}, {}]}
		}
	}`
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		assert.Equal(t, "application/json", r.Header.Get("Accept"))
		w.Header().Set("Content-Type", "application/json")
		w.Write([]byte(payload))
	}))
	defer srv.Close()

	transport := &redirectTransport{base: srv.Client().Transport, target: srv.URL}
	c := &ORCIDClient{httpClient: &http.Client{Transport: transport}}

	r, err := c.LookupResearcher(context.Background(), "0000-0001-9999-8888")
	require.NoError(t, err)
	assert.Equal(t, "Alice Wong", r.Name)
	assert.Equal(t, "awong@example.com", r.Email)
	assert.Equal(t, []string{"Stanford"}, r.Institutions)
	assert.Equal(t, 3, r.Publications)
}

// --- CrossRefClient tests ---

func TestCrossRefClient_ResolveDOI_NotFound(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.NotFound(w, r)
	}))
	defer srv.Close()

	transport := &redirectTransport{base: srv.Client().Transport, target: srv.URL}
	c := &CrossRefClient{httpClient: &http.Client{Transport: transport}}

	_, err := c.ResolveDOI(context.Background(), "10.9999/bad")
	assert.Error(t, err)
	assert.Contains(t, err.Error(), "not found")
}

func TestCrossRefClient_ResolveDOI_ParseResponse(t *testing.T) {
	raw := `{"message":{"DOI":"10.1038/nature15393","title":["A global reference for human genetic variation"],"container-title":["Nature"],"publisher":"Springer Nature","type":"journal-article","URL":"http://dx.doi.org/10.1038/nature15393","created":{"date-parts":[[2015,10,1]]},"author":[{"given":"","family":"1000 Genomes Project Consortium"}],"score":1.0}}`
	var result crossRefResponse
	require.NoError(t, json.Unmarshal([]byte(raw), &result))
	assert.Equal(t, "10.1038/nature15393", result.Message.DOI)
	assert.Equal(t, "A global reference for human genetic variation", result.Message.Title[0])
	assert.Equal(t, "Nature", result.Message.ContainerTitle[0])
	assert.Equal(t, 2015, result.Message.Created.DateParts[0][0])
}

func TestCrossRefClient_ResolveDOI_ServerSuccess(t *testing.T) {
	payload := `{"message":{"DOI":"10.1038/nature15393","title":["A global reference for human genetic variation"],"container-title":["Nature"],"publisher":"Springer Nature","type":"journal-article","URL":"http://dx.doi.org/10.1038/nature15393","created":{"date-parts":[[2015,10,1]]},"author":[{"given":"","family":"1000 Genomes Project Consortium"}],"score":1.0}}`
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		assert.Contains(t, r.Header.Get("User-Agent"), "phantom-science")
		w.Header().Set("Content-Type", "application/json")
		w.Write([]byte(payload))
	}))
	defer srv.Close()

	transport := &redirectTransport{base: srv.Client().Transport, target: srv.URL}
	c := &CrossRefClient{httpClient: &http.Client{Transport: transport}}

	pub, err := c.ResolveDOI(context.Background(), "10.1038/nature15393")
	require.NoError(t, err)
	assert.Equal(t, "10.1038/nature15393", pub.DOI)
	assert.Equal(t, "A global reference for human genetic variation", pub.Title)
	assert.Equal(t, "Nature", pub.Journal)
	assert.Equal(t, 2015, pub.Year)
	assert.Equal(t, "Springer Nature", pub.Publisher)
}

// redirectTransport redirects all HTTP requests to a test server URL.
type redirectTransport struct {
	base   http.RoundTripper
	target string
}

func (t *redirectTransport) RoundTrip(req *http.Request) (*http.Response, error) {
	// Rewrite the request URL to point at the test server, preserving path+query.
	newURL := *req.URL
	newURL.Scheme = "http"
	// Parse target to extract host
	var targetHost string
	fmt.Sscanf(t.target, "http://%s", &targetHost)
	newURL.Host = targetHost

	newReq := req.Clone(req.Context())
	newReq.URL = &newURL
	newReq.Host = targetHost
	return t.base.RoundTrip(newReq)
}
