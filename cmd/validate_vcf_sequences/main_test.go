//go:build !release
// +build !release

package main

import "testing"

func TestCreateVCFValidator_ReturnsNil(t *testing.T) {
	err := createVCFValidator()
	if err != nil {
		t.Errorf("createVCFValidator() = %v, want nil", err)
	}
}
