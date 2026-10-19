//go:build !release

package reference_test

import (
	"fmt"
	"os"
	"sync"
	"testing"
	"time"

	"github.com/scttfrdmn/phantom-science/internal/genomics/reference"
	"github.com/stretchr/testify/require"
)

// TestReferenceManager_OptimizeCache_ConcurrentWriter: OptimizeCache must not
// take loaderMutex.RLock twice (it calls GetCombinedStatistics, which also
// read-locks). A second RLock blocks behind a pending Lock from
// RegisterReferenceGenome, which in turn waits for the first RLock: a
// deadlock. Run with -race to also check the request counters.
func TestReferenceManager_OptimizeCache_ConcurrentWriter(t *testing.T) {
	p := makeDummyFastaFile(t)
	rm, err := reference.NewReferenceManager(reference.ManagerConfig{
		DefaultReferenceFile: p,
		DefaultCacheSize:     1 << 20,
		MaxLoaders:           1 << 30,
	})
	require.NoError(t, err)

	// OptimizeCache prints suggestions; keep the test output readable.
	devNull, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	require.NoError(t, err)
	stdout := os.Stdout
	os.Stdout = devNull
	defer func() { os.Stdout = stdout; devNull.Close() }()

	const iters = 2000
	done := make(chan struct{})
	go func() {
		defer close(done)
		var wg sync.WaitGroup
		wg.Add(4)
		go func() {
			defer wg.Done()
			for i := 0; i < iters; i++ {
				rm.OptimizeCache()
			}
		}()
		go func() {
			defer wg.Done()
			for i := 0; i < iters; i++ {
				_ = rm.RegisterReferenceGenome(fmt.Sprintf("ref%d", i), p, 1024)
			}
		}()
		go func() {
			defer wg.Done()
			for i := 0; i < iters; i++ {
				_, _ = rm.GetSequence("chr1", 0, 4)
			}
		}()
		go func() {
			defer wg.Done()
			for i := 0; i < iters; i++ {
				_ = rm.GetCombinedStatistics()
			}
		}()
		wg.Wait()
	}()

	select {
	case <-done:
	case <-time.After(30 * time.Second):
		os.Stdout = stdout
		t.Fatal("OptimizeCache deadlocked against a concurrent RegisterReferenceGenome")
	}
}
