// Package blast provides BLAST+ integration for sequence alignment.
package blast

import (
	"fmt"
	"os/exec"
)

// BinaryPath locates a BLAST+ executable in the system PATH.
// It returns an error if the tool is not installed.
func BinaryPath(name string) (string, error) {
	path, err := exec.LookPath(name)
	if err != nil {
		return "", fmt.Errorf("BLAST+ tool '%s' not found in PATH. Please install NCBI BLAST+ (including %s) or run 'make download-blast' to fetch it.", name, name)
	}
	return path, nil
}
