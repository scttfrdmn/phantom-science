//go:build !release

package s3sink

import (
	"bytes"
	"context"
	"io"
	"math/rand/v2"
	"net/url"
	"os"
	"path/filepath"
	"sort"
	"sync"
	"testing"

	"github.com/aws/aws-sdk-go/aws"
	"github.com/aws/aws-sdk-go/aws/credentials"
	"github.com/aws/aws-sdk-go/aws/request"
	"github.com/aws/aws-sdk-go/aws/session"
	"github.com/aws/aws-sdk-go/service/s3"
	"github.com/aws/aws-sdk-go/service/s3/s3manager"
)

func TestParseURI(t *testing.T) {
	for in, want := range map[string]Location{
		"s3://b":             {Bucket: "b"},
		"s3://b/":            {Bucket: "b"},
		"s3://b/run1":        {Bucket: "b", Prefix: "run1/"},
		"s3://b/a/b/c/":      {Bucket: "b", Prefix: "a/b/c/"},
		"s3://my.bucket/x/y": {Bucket: "my.bucket", Prefix: "x/y/"},
	} {
		got, err := ParseURI(in)
		if err != nil || got != want {
			t.Errorf("ParseURI(%q) = %+v, %v; want %+v", in, got, err, want)
		}
	}
	for _, bad := range []string{"", "/tmp/out", "s3://", "s3:///prefix", "s3://b/a//c", "s3://b/../x", "s3://b/./x"} {
		if _, err := ParseURI(bad); err == nil {
			t.Errorf("ParseURI(%q) accepted", bad)
		}
	}
	l, _ := ParseURI("s3://b/run")
	if got := l.Key("FC1/x.fastq.gz"); got != "run/FC1/x.fastq.gz" {
		t.Errorf("Key = %q", got)
	}
	if got := l.URI(); got != "s3://b/run/" {
		t.Errorf("URI = %q", got)
	}
}

func TestEncodeTags(t *testing.T) {
	got := EncodeTags(map[string]string{"b": "2 3", "a": "x&y"})
	if got != "a=x%26y&b=2+3" {
		t.Errorf("EncodeTags = %q", got)
	}
	v, err := url.ParseQuery(got)
	if err != nil || v.Get("a") != "x&y" || v.Get("b") != "2 3" {
		t.Errorf("tags do not round-trip: %v %v", v, err)
	}
}

type recordingUploader struct {
	keys []string
	data map[string][]byte
	err  error
}

func (r *recordingUploader) Upload(_ context.Context, in UploadInput) error {
	b, err := io.ReadAll(in.Body)
	if err != nil {
		return err
	}
	if int64(len(b)) != in.Size {
		return io.ErrShortWrite
	}
	if r.data == nil {
		r.data = map[string][]byte{}
	}
	r.keys = append(r.keys, in.Bucket+"/"+in.Key)
	r.data[in.Key] = b
	return r.err
}

func TestSink_PutUploadsAndRemovesSpool(t *testing.T) {
	dir := t.TempDir()
	local := filepath.Join(dir, "a.fastq")
	if err := os.WriteFile(local, []byte("@r\nACGT\n+\nIIII\n"), 0o600); err != nil {
		t.Fatal(err)
	}
	up := &recordingUploader{}
	s := New(Location{Bucket: "bk", Prefix: "p/"}, up)
	key, err := s.Put(context.Background(), local, "FC/a.fastq", nil)
	if err != nil || key != "p/FC/a.fastq" {
		t.Fatalf("Put = %q, %v", key, err)
	}
	if len(up.keys) != 1 || up.keys[0] != "bk/p/FC/a.fastq" || string(up.data[key]) != "@r\nACGT\n+\nIIII\n" {
		t.Errorf("uploaded %v", up.keys)
	}
	if _, err := os.Stat(local); !os.IsNotExist(err) {
		t.Error("spooled file not removed after upload")
	}

	// A failed upload is reported and the spool is still cleaned up.
	if err := os.WriteFile(local, []byte("x"), 0o600); err != nil {
		t.Fatal(err)
	}
	up.err = io.ErrUnexpectedEOF
	if _, err := s.Put(context.Background(), local, "a.fastq", nil); err == nil {
		t.Error("upload error not returned")
	}
	if _, err := os.Stat(local); !os.IsNotExist(err) {
		t.Error("spooled file not removed after a failed upload")
	}
}

// fakeS3 implements the multipart calls s3manager makes; the embedded
// client (never used for network calls) supplies request presigning.
type fakeS3 struct {
	*s3.S3
	mu      sync.Mutex
	created []*s3.CreateMultipartUploadInput
	parts   map[int64][]byte
	done    bool
	puts    int
}

func (f *fakeS3) CreateMultipartUploadWithContext(_ aws.Context, in *s3.CreateMultipartUploadInput, _ ...request.Option) (*s3.CreateMultipartUploadOutput, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	f.created = append(f.created, in)
	return &s3.CreateMultipartUploadOutput{UploadId: aws.String("up-1")}, nil
}

func (f *fakeS3) UploadPartWithContext(_ aws.Context, in *s3.UploadPartInput, _ ...request.Option) (*s3.UploadPartOutput, error) {
	b, err := io.ReadAll(in.Body)
	if err != nil {
		return nil, err
	}
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.parts == nil {
		f.parts = map[int64][]byte{}
	}
	f.parts[*in.PartNumber] = b
	return &s3.UploadPartOutput{ETag: aws.String("etag")}, nil
}

func (f *fakeS3) CompleteMultipartUploadWithContext(_ aws.Context, in *s3.CompleteMultipartUploadInput, _ ...request.Option) (*s3.CompleteMultipartUploadOutput, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	f.done = len(in.MultipartUpload.Parts) == len(f.parts)
	return &s3.CompleteMultipartUploadOutput{}, nil
}

func (f *fakeS3) PutObjectRequest(in *s3.PutObjectInput) (*request.Request, *s3.PutObjectOutput) {
	f.mu.Lock()
	f.puts++
	f.mu.Unlock()
	return f.S3.PutObjectRequest(in)
}

// The real uploader sends a file larger than one part as a multipart
// upload whose parts reassemble the file, with the object tags.
func TestNewUploader_Multipart(t *testing.T) {
	sess, err := session.NewSession(&aws.Config{
		Region:      aws.String("us-east-1"),
		Credentials: credentials.NewStaticCredentials("id", "secret", ""),
		Endpoint:    aws.String("http://127.0.0.1:1"),
	})
	if err != nil {
		t.Fatal(err)
	}
	fake := &fakeS3{S3: s3.New(sess)}
	up := NewUploader(fake, s3manager.MinUploadPartSize)

	content := make([]byte, 2*s3manager.MinUploadPartSize+12345)
	r := rand.New(rand.NewPCG(1, 2))
	for i := range content {
		content[i] = "ACGT"[r.IntN(4)]
	}
	local := filepath.Join(t.TempDir(), "big.fastq")
	if err := os.WriteFile(local, content, 0o600); err != nil {
		t.Fatal(err)
	}
	s := New(Location{Bucket: "bk", Prefix: "run/"}, up)
	key, err := s.Put(context.Background(), local, "big.fastq", map[string]string{"phantom:seed": "7"})
	if err != nil {
		t.Fatal(err)
	}
	if fake.puts != 0 || len(fake.created) != 1 || !fake.done {
		t.Fatalf("puts=%d creates=%d completed=%v: want one completed multipart upload", fake.puts, len(fake.created), fake.done)
	}
	c := fake.created[0]
	if *c.Bucket != "bk" || *c.Key != key || key != "run/big.fastq" {
		t.Errorf("multipart upload to %s/%s, key %s", *c.Bucket, *c.Key, key)
	}
	if c.Tagging == nil || *c.Tagging != "phantom%3Aseed=7" {
		t.Errorf("tagging = %v", c.Tagging)
	}
	nums := make([]int64, 0, len(fake.parts))
	for n := range fake.parts {
		nums = append(nums, n)
	}
	sort.Slice(nums, func(i, j int) bool { return nums[i] < nums[j] })
	var joined []byte
	for _, n := range nums {
		joined = append(joined, fake.parts[n]...)
	}
	if len(nums) < 3 || !bytes.Equal(joined, content) {
		t.Errorf("%d parts reassemble to %d bytes, want %d bytes in >= 3 parts", len(nums), len(joined), len(content))
	}
	if _, err := os.Stat(local); !os.IsNotExist(err) {
		t.Error("spooled file not removed")
	}
}
