// Package fasta provides FASTA file reading and writing utilities.
package fasta

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

// Record represents a single FASTA record.
type Record struct {
	ID   string // The sequence identifier (first word of header)
	Desc string // The full description (everything after the ID in header)
	Seq  string // The sequence string (uppercase)
}

// ReadAll reads all FASTA records from an io.Reader.
func ReadAll(r io.Reader) ([]Record, error) {
	var records []Record
	var rec *Record
	var seqBuilder strings.Builder

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024) // support long lines

	for scanner.Scan() {
		line := strings.TrimRight(scanner.Text(), "\r\n")
		if strings.HasPrefix(line, ">") {
			if rec != nil {
				rec.Seq = strings.ToUpper(seqBuilder.String())
				records = append(records, *rec)
				seqBuilder.Reset()
			}
			header := line[1:]
			parts := strings.SplitN(header, " ", 2)
			rec = &Record{}
			rec.ID = parts[0]
			if len(parts) > 1 {
				rec.Desc = parts[1]
			}
		} else if rec != nil {
			seqBuilder.WriteString(line)
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	if rec != nil {
		rec.Seq = strings.ToUpper(seqBuilder.String())
		records = append(records, *rec)
	}
	return records, nil
}

// ReadFile reads all FASTA records from a file.
func ReadFile(path string) ([]Record, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ReadAll(f)
}

// Writer writes FASTA records to an io.Writer.
type Writer struct {
	w       io.Writer
	lineLen int
}

// NewWriter creates a new FASTA writer with default 60 chars per line.
func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w, lineLen: 60}
}

// NewWriterWidth creates a new FASTA writer with specified line width.
func NewWriterWidth(w io.Writer, width int) *Writer {
	return &Writer{w: w, lineLen: width}
}

// Write writes a single FASTA record.
func (w *Writer) Write(rec *Record) error {
	header := rec.ID
	if rec.Desc != "" {
		header += " " + rec.Desc
	}
	if _, err := fmt.Fprintf(w.w, ">%s\n", header); err != nil {
		return err
	}
	seq := rec.Seq
	for i := 0; i < len(seq); i += w.lineLen {
		end := i + w.lineLen
		if end > len(seq) {
			end = len(seq)
		}
		if _, err := fmt.Fprintf(w.w, "%s\n", seq[i:end]); err != nil {
			return err
		}
	}
	return nil
}

// WriteAll writes all FASTA records.
func (w *Writer) WriteAll(records []Record) error {
	for i := range records {
		if err := w.Write(&records[i]); err != nil {
			return err
		}
	}
	return nil
}

// WriteFile writes all FASTA records to a file.
func WriteFile(path string, records []Record) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := NewWriter(f)
	return w.WriteAll(records)
}
