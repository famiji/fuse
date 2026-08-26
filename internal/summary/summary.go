// Package summary provides summary table generation for Fuse.
// It combines multiple abricate report files into a presence/absence matrix.
package summary

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strings"
)

// Absent is the symbol for absent genes.
const Absent = "."

// FieldSep is the separator between multiple values in a cell.
const FieldSep = ";"

// Options controls summary table formatting.
type Options struct {
	UseIdentity bool // Use %IDENTITY instead of %COVERAGE
	NoPath      bool
	CSV         bool
}

// Generate reads multiple report files and produces a summary table.
func Generate(filenames []string, opts Options, w *bufio.Writer) error {
	sep := "\t"
	if opts.CSV {
		sep = ","
	}

	field := "%COVERAGE"
	if opts.UseIdentity {
		field = "%IDENTITY"
	}

	// Dutch mode: single file, use #FILE column as key
	dutchMode := len(filenames) == 1

	geneSet := make(map[string]bool)
	data := make(map[string]map[string][]map[string]string)
	var header []string

	for _, fname := range filenames {
		if _, exists := data[fname]; exists {
			fmt.Fprintf(os.Stderr, "WARNING: Skipping duplicate file: %s\n", fname)
			continue
		}
		if !dutchMode {
			data[fname] = make(map[string][]map[string]string)
		}

		f, err := os.Open(fname)
		if err != nil {
			return fmt.Errorf("can't open '%s': %w", fname, err)
		}

		scanner := bufio.NewScanner(f)
		scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
		for scanner.Scan() {
			line := scanner.Text()
			cols := strings.Split(line, sep)
			if len(cols) == 0 {
				continue
			}
			if header == nil {
				header = cols
			}
			if strings.HasPrefix(cols[0], "#") {
				continue
			}

			file := fname
			if dutchMode {
				file = cols[0]
			}
			if opts.NoPath {
				file = filepath.Base(file)
				if dutchMode {
					cols[0] = filepath.Base(cols[0])
					file = cols[0]
				}
			}

			row := make(map[string]string)
			for i, val := range cols {
				if i < len(header) {
					row[header[i]] = val
				}
			}

			gene := row["GENE"]
			geneSet[gene] = true
			if data[file] == nil {
				data[file] = make(map[string][]map[string]string)
			}
			data[file][gene] = append(data[file][gene], row)
		}
		f.Close()
		if err := scanner.Err(); err != nil {
			return err
		}
	}

	// Sort genes
	var genes []string
	for g := range geneSet {
		genes = append(genes, g)
	}
	slices.Sort(genes)

	// Print header
	outHeader := []string{"#FILE", "NUM_FOUND"}
	outHeader = append(outHeader, genes...)
	fmt.Fprintln(w, strings.Join(outHeader, sep))

	// Sort file names
	var fileNames []string
	for f := range data {
		fileNames = append(fileNames, f)
	}
	slices.Sort(fileNames)

	for _, fname := range fileNames {
		var present []string
		for _, gene := range genes {
			if rows, exists := data[fname][gene]; exists {
				var values []string
				for _, r := range rows {
					values = append(values, r[field])
				}
				present = append(present, strings.Join(values, FieldSep))
			} else {
				present = append(present, Absent)
			}
		}
		label := fname
		if opts.NoPath {
			label = filepath.Base(label)
		}
		numFound := len(data[fname])
		row := []string{label, fmt.Sprintf("%d", numFound)}
		row = append(row, present...)
		fmt.Fprintln(w, strings.Join(row, sep))
	}

	return nil
}
