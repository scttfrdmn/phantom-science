// Package s3sink uploads finished output files to S3 (#16).
//
// The engine writes each output file to a local spool directory, finalizes
// it there (fingerprint embedding, compression), then hands it to a Sink,
// which uploads it under the run's prefix (multipart for large files) and
// deletes the spooled copy. Keys mirror the local file names: an output
// that would be written as <dir>/<rel> is uploaded as <prefix><rel>.
//
// Uploads go through the small Uploader interface so tests use a fake; the
// real implementation (NewUploader, NewDefaultUploader) is aws-sdk-go v1's
// s3manager, with credentials and region from the default chain.
//
// Testing is limited to fakes: the Sink tests use a fake Uploader, and the
// s3manager path is exercised against a fake S3 client that records the
// multipart calls. Nothing in the repository runs against real AWS (or an
// S3-compatible store), so permissions, tagging limits, retries and other
// service-side behaviour are untested.
package s3sink

import (
	"context"
	"errors"
	"fmt"
	"net/url"
	"os"
	"path"
	"sort"
	"strings"

	"github.com/aws/aws-sdk-go/aws"
	"github.com/aws/aws-sdk-go/aws/session"
	"github.com/aws/aws-sdk-go/service/s3"
	"github.com/aws/aws-sdk-go/service/s3/s3iface"
	"github.com/aws/aws-sdk-go/service/s3/s3manager"
)

// Scheme is the URI scheme of S3 destinations.
const Scheme = "s3://"

// DefaultPartSize is the multipart part size of the real uploader. S3
// allows 10,000 parts, so this supports objects up to about 640 GB.
const DefaultPartSize = 64 << 20

// Location is an S3 bucket and key prefix.
type Location struct {
	Bucket string
	// Prefix is empty or ends in "/".
	Prefix string
}

// IsURI reports whether s is an s3:// URI.
func IsURI(s string) bool { return strings.HasPrefix(s, Scheme) }

// ParseURI parses s3://bucket[/prefix[/]]. The prefix is treated as a
// directory: a trailing "/" is added when missing.
func ParseURI(s string) (Location, error) {
	if !IsURI(s) {
		return Location{}, fmt.Errorf("%q is not an s3:// URI", s)
	}
	bucket, prefix, _ := strings.Cut(strings.TrimPrefix(s, Scheme), "/")
	if bucket == "" {
		return Location{}, fmt.Errorf("S3 URI %q has no bucket", s)
	}
	prefix = strings.TrimLeft(prefix, "/")
	if prefix != "" && !strings.HasSuffix(prefix, "/") {
		prefix += "/"
	}
	for _, seg := range strings.Split(strings.TrimSuffix(prefix, "/"), "/") {
		if seg == "." || seg == ".." || (seg == "" && prefix != "") {
			return Location{}, fmt.Errorf("S3 URI %q: invalid prefix", s)
		}
	}
	return Location{Bucket: bucket, Prefix: prefix}, nil
}

// URI returns the location as s3://bucket/prefix.
func (l Location) URI() string { return Scheme + l.Bucket + "/" + l.Prefix }

// Key returns the object key of rel (a slash-separated relative path).
func (l Location) Key(rel string) string { return l.Prefix + path.Clean(rel) }

// UploadInput is one object upload.
type UploadInput struct {
	Bucket, Key string
	// Body is the object content; Size its length in bytes.
	Body *os.File
	Size int64
	// Tags are set as S3 object tags.
	Tags map[string]string
}

// Uploader uploads one object.
type Uploader interface {
	Upload(ctx context.Context, in UploadInput) error
}

// Sink uploads spooled files under a Location.
type Sink struct {
	loc Location
	up  Uploader
}

// New returns a Sink uploading to loc with up.
func New(loc Location, up Uploader) *Sink { return &Sink{loc: loc, up: up} }

// Location returns the sink's destination.
func (s *Sink) Location() Location { return s.loc }

// Put uploads the local file to the key of rel and removes the local file
// (also when the upload fails: the spool is temporary). It returns the key.
func (s *Sink) Put(ctx context.Context, local, rel string, tags map[string]string) (string, error) {
	key := s.loc.Key(rel)
	f, err := os.Open(local) // #nosec G304 -- a spooled output file of this run
	if err != nil {
		return "", err
	}
	info, err := f.Stat()
	if err == nil {
		err = s.up.Upload(ctx, UploadInput{Bucket: s.loc.Bucket, Key: key, Body: f, Size: info.Size(), Tags: tags})
		if err != nil {
			err = fmt.Errorf("uploading s3://%s/%s: %w", s.loc.Bucket, key, err)
		}
	}
	err = errors.Join(err, f.Close(), os.Remove(local))
	if err != nil {
		return "", err
	}
	return key, nil
}

// EncodeTags returns tags as an S3 tagging query string, keys sorted.
func EncodeTags(tags map[string]string) string {
	keys := make([]string, 0, len(tags))
	for k := range tags {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	var b strings.Builder
	for i, k := range keys {
		if i > 0 {
			b.WriteByte('&')
		}
		b.WriteString(url.QueryEscape(k))
		b.WriteByte('=')
		b.WriteString(url.QueryEscape(tags[k]))
	}
	return b.String()
}

// managerUploader is the s3manager-backed Uploader.
type managerUploader struct{ u *s3manager.Uploader }

// NewUploader returns an Uploader that uses s3manager over client:
// multipart uploads of partSize bytes (DefaultPartSize when <= 0) for
// objects larger than one part.
func NewUploader(client s3iface.S3API, partSize int64) Uploader {
	if partSize <= 0 {
		partSize = DefaultPartSize
	}
	return managerUploader{s3manager.NewUploaderWithClient(client, func(u *s3manager.Uploader) {
		u.PartSize = max(partSize, s3manager.MinUploadPartSize)
	})}
}

// NewDefaultUploader returns the s3manager Uploader with the default AWS
// credential chain and region (environment, shared config, instance role).
func NewDefaultUploader() (Uploader, error) {
	sess, err := session.NewSessionWithOptions(session.Options{SharedConfigState: session.SharedConfigEnable})
	if err != nil {
		return nil, fmt.Errorf("creating AWS session: %w", err)
	}
	return NewUploader(s3.New(sess), DefaultPartSize), nil
}

// Upload implements Uploader.
func (m managerUploader) Upload(ctx context.Context, in UploadInput) error {
	input := &s3manager.UploadInput{Bucket: aws.String(in.Bucket), Key: aws.String(in.Key), Body: in.Body}
	if len(in.Tags) > 0 {
		input.Tagging = aws.String(EncodeTags(in.Tags))
	}
	_, err := m.u.UploadWithContext(ctx, input)
	return err
}
