// Package seqio provides sequence file I/O utilities including
// decompression of gzip/bzip2 files and conversion from GenBank/EMBL to FASTA.
package seqio

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/fmj/fuse/internal/fasta"
)

// AnyToFASTA converts any supported input file to FASTA format.
// It handles gzip and bzip2 compression, and GenBank/EMBL format conversion.
// If any2fasta is available on PATH, it uses that; otherwise it falls back
// to built-in conversion.
func AnyToFASTA(path string) ([]fasta.Record, error) {
	r, cleanup, err := openCompressed(path)
	if err != nil {
		return nil, err
	}
	defer cleanup()

	// Check if it's already FASTA
	br := bufio.NewReader(r)
	firstLine, err := br.ReadString('\n')
	if err != nil && err != io.EOF {
		return nil, err
	}
	firstLine = strings.TrimRight(firstLine, "\r\n")

	if strings.HasPrefix(firstLine, ">") {
		// It's FASTA - read the rest
		var buf bytes.Buffer
		buf.WriteString(firstLine)
		buf.WriteString("\n")
		_, err = io.Copy(&buf, br)
		if err != nil {
			return nil, err
		}
		return fasta.ReadAll(&buf)
	}

	// Check for GenBank or EMBL
	upperFirst := strings.ToUpper(strings.TrimSpace(firstLine))
	if strings.HasPrefix(upperFirst, "LOCUS") {
		return parseGenBank(br, firstLine)
	}
	if strings.HasPrefix(upperFirst, "ID") || strings.HasPrefix(upperFirst, "ID ") {
		return parseEMBL(br, firstLine)
	}

	// Try any2fasta if available
	if exe, _ := exec.LookPath("any2fasta"); exe != "" {
		return convertWithAny2Fasta(path)
	}

	return nil, fmt.Errorf("unsupported file format (first line: %q)", firstLine)
}

// openCompressed opens a file, transparently decompressing gzip or bzip2.
func openCompressed(path string) (io.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}

	ext := strings.ToLower(filepath.Ext(path))
	switch ext {
	case ".gz":
		gr, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, nil, err
		}
		cleanup := func() {
			gr.Close()
			f.Close()
		}
		return gr, cleanup, nil
	case ".bz2":
		// Use bzip2 command since Go stdlib doesn't have bzip2 decompressor
		// Actually, we can use the bzip2 decompression via exec
		cmd := exec.Command("bzip2", "-dc", path)
		stdout, err := cmd.StdoutPipe()
		if err != nil {
			f.Close()
			return nil, nil, err
		}
		if err := cmd.Start(); err != nil {
			f.Close()
			return nil, nil, err
		}
		cleanup := func() {
			stdout.Close()
			cmd.Wait()
			f.Close()
		}
		return stdout, cleanup, nil
	default:
		cleanup := func() { f.Close() }
		return f, cleanup, nil
	}
}

// convertWithAny2Fasta uses the external any2fasta tool.
func convertWithAny2Fasta(path string) ([]fasta.Record, error) {
	cmd := exec.Command("any2fasta", "-q", "-u", path)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		return nil, fmt.Errorf("any2fasta failed: %w: %s", err, stderr.String())
	}
	return fasta.ReadAll(&stdout)
}

// parseGenBank parses a GenBank file and extracts FASTA sequences.
func parseGenBank(r io.Reader, firstLine string) ([]fasta.Record, error) {
	var records []fasta.Record
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)

	// We need to re-read from the beginning including firstLine
	var content strings.Builder
	content.WriteString(firstLine)
	content.WriteString("\n")
	for scanner.Scan() {
		content.WriteString(scanner.Text())
		content.WriteString("\n")
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	lines := strings.Split(content.String(), "\n")

	var currentID string
	var currentDesc string
	var seqBuilder strings.Builder
	inOrigin := false

	flush := func() {
		if currentID != "" && seqBuilder.Len() > 0 {
			records = append(records, fasta.Record{
				ID:   currentID,
				Desc: currentDesc,
				Seq:  strings.ToUpper(strings.ReplaceAll(seqBuilder.String(), " ", "")),
			})
		}
		currentID = ""
		currentDesc = ""
		seqBuilder.Reset()
	}

	for _, line := range lines {
		if strings.HasPrefix(line, "LOCUS") {
			flush()
			inOrigin = false
			// Parse LOCUS line: LOCUS   LGJG01000042  384079 bp  DNA  linear  BCT 24-SEP-2015
			fields := strings.Fields(line)
			if len(fields) >= 2 {
				currentID = fields[1]
			}
			// Build description from rest
			if len(fields) > 2 {
				currentDesc = strings.Join(fields[2:], " ")
			}
		} else if strings.HasPrefix(line, "DEFINITION") {
			// DEFINITION  Staphylococcus epidermidis strain ST14 contig_1...
			desc := strings.TrimSpace(strings.TrimPrefix(line, "DEFINITION"))
			desc = strings.TrimRight(desc, ".")
			if currentDesc == "" {
				currentDesc = desc
			} else {
				currentDesc = desc
			}
		} else if strings.HasPrefix(line, "ORIGIN") {
			inOrigin = true
		} else if line == "//" {
			flush()
			inOrigin = false
		} else if inOrigin {
			// Sequence lines: "  1 aagctgact..."
			// Remove leading numbers and spaces
			trimmed := strings.TrimSpace(line)
			// Remove leading position number
			parts := strings.Fields(trimmed)
			for _, p := range parts {
				if p != "" && !isAllDigits(p) {
					seqBuilder.WriteString(p)
				}
			}
		}
	}
	flush()

	if len(records) == 0 {
		return nil, fmt.Errorf("no sequences found in GenBank file")
	}
	return records, nil
}

func isAllDigits(s string) bool {
	for _, c := range s {
		if c < '0' || c > '9' {
			return false
		}
	}
	return true
}

// parseEMBL parses an EMBL file and extracts FASTA sequences.
func parseEMBL(r io.Reader, firstLine string) ([]fasta.Record, error) {
	var records []fasta.Record
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)

	var content strings.Builder
	content.WriteString(firstLine)
	content.WriteString("\n")
	for scanner.Scan() {
		content.WriteString(scanner.Text())
		content.WriteString("\n")
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	lines := strings.Split(content.String(), "\n")

	var currentID string
	var currentDesc string
	var seqBuilder strings.Builder
	inSQ := false

	flush := func() {
		if currentID != "" && seqBuilder.Len() > 0 {
			records = append(records, fasta.Record{
				ID:   currentID,
				Desc: currentDesc,
				Seq:  strings.ToUpper(strings.ReplaceAll(seqBuilder.String(), " ", "")),
			})
		}
		currentID = ""
		currentDesc = ""
		seqBuilder.Reset()
	}

	for _, line := range lines {
		if strings.HasPrefix(line, "ID") && !strings.HasPrefix(line, "ID ") == false {
			flush()
			inSQ = false
			// ID   X12345; SV 1; linear; genomic DNA; STD; PRO; 384079 BP.
			fields := strings.Fields(line)
			if len(fields) >= 2 {
				currentID = strings.TrimRight(fields[1], ";")
			}
		} else if strings.HasPrefix(line, "DE") {
			desc := strings.TrimSpace(strings.TrimPrefix(line, "DE"))
			desc = strings.TrimPrefix(desc, "DE")
			desc = strings.TrimSpace(desc)
			if currentDesc == "" {
				currentDesc = desc
			}
		} else if strings.HasPrefix(line, "SQ") {
			inSQ = true
		} else if line == "//" {
			flush()
			inSQ = false
		} else if inSQ {
			trimmed := strings.TrimSpace(line)
			parts := strings.Fields(trimmed)
			for _, p := range parts {
				if p != "" && !isAllDigits(p) {
					seqBuilder.WriteString(p)
				}
			}
		}
	}
	flush()

	if len(records) == 0 {
		return nil, fmt.Errorf("no sequences found in EMBL file")
	}
	return records, nil
}

// IsCompressed returns true if the file appears to be compressed.
func IsCompressed(path string) bool {
	ext := strings.ToLower(filepath.Ext(path))
	return ext == ".gz" || ext == ".bz2"
}
