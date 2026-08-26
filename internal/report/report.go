// Package report provides result formatting and output for Fuse.
package report

import (
	"cmp"
	"fmt"
	"path/filepath"
	"slices"
	"strconv"
	"strings"

	"github.com/fmj/fuse/internal/blast"
	"github.com/fmj/fuse/internal/database"
)

// Columns are the output column headers.
var Columns = []string{
	"#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE",
	"COVERAGE", "COVERAGE_MAP", "GAPS",
	"%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION",
	"PRODUCT", "RESISTANCE",
}

// Result represents a single hit in the output format.
type Result struct {
	File        string
	Sequence    string
	Start       int
	End         int
	Strand      string
	Gene        string
	Coverage    string
	CoverageMap string
	Gaps        string
	PcCoverage  float64
	PcIdentity  float64
	Database    string
	Accession   string
	Product     string
	Resistance  string
}

// Options controls report formatting.
type Options struct {
	NoPath   bool
	CSV      bool
	NoHeader bool
}

// FormatResult converts a BLAST hit into a Result row.
func FormatResult(hit blast.Hit, fileName string, dbName string, minCov float64, noPath bool) (*Result, bool) {
	// Calculate coverage percentage
	pccov := 100.0 * float64(hit.Length-hit.Gaps) / float64(hit.SLen)
	if pccov < minCov {
		return nil, false
	}

	// Parse the subject sequence ID
	db, gene, acc, resistance := database.ParseSeqID(hit.SSeqID, dbName)
	if gene == "" && acc == "" {
		acc = ""
		gene = hit.SSeqID
		db = dbName
	}

	// Product
	product := hit.STitle
	if product == "" {
		product = "n/a"
	}
	// Remove commas and tabs from product
	product = strings.ReplaceAll(product, ",", "")
	product = strings.ReplaceAll(product, "\t", "")
	// If product contains IDSeparator, strip the first word
	if strings.Contains(product, database.IDSeparator) {
		parts := strings.SplitN(product, " ", 2)
		if len(parts) > 1 {
			product = parts[1]
		}
	}

	// Coverage map
	miniMap := MiniMap(hit.SStart, hit.SEnd, hit.SLen, hit.GapOpen)

	// File path
	file := fileName
	if noPath {
		file = filepath.Base(file)
	}

	return &Result{
		File:        file,
		Sequence:    hit.QSeqID,
		Start:       hit.QStart,
		End:         hit.QEnd,
		Strand:      strandSymbol(hit.SStrand),
		Gene:        gene,
		Coverage:    fmt.Sprintf("%d-%d/%d", hit.SStart, hit.SEnd, hit.SLen),
		CoverageMap: miniMap,
		Gaps:        fmt.Sprintf("%d/%d", hit.GapOpen, hit.Gaps),
		PcCoverage:  pccov,
		PcIdentity:  hit.PIdent,
		Database:    db,
		Accession:   acc,
		Product:     product,
		Resistance:  resistance,
	}, true
}

func strandSymbol(s string) string {
	if s == "minus" {
		return "-"
	}
	return "+"
}

// MiniMap creates a visual representation of the hit coverage.
func MiniMap(x, y, l, broken int) string {
	width := 15
	if broken > 0 {
		width = 14
	}
	scale := float64(l) / float64(width)
	if scale == 0 {
		scale = 1
	}
	xScaled := int(float64(x) / scale)
	yScaled := int(float64(y) / scale)
	lScaled := int(float64(l) / scale)

	var sb strings.Builder
	for i := 0; i < width; i++ {
		if i >= xScaled && i <= yScaled {
			sb.WriteByte('=')
		} else {
			sb.WriteByte('.')
		}
		if broken > 0 && i == width/2 {
			sb.WriteByte('/')
		}
	}
	_ = lScaled // suppress unused warning
	return sb.String()
}

// FormatLine formats a Result as a tab or comma separated line.
func FormatLine(r *Result, sep string) string {
	fields := []string{
		r.File,
		r.Sequence,
		strconv.Itoa(r.Start),
		strconv.Itoa(r.End),
		r.Strand,
		r.Gene,
		r.Coverage,
		r.CoverageMap,
		r.Gaps,
		fmt.Sprintf("%.2f", r.PcCoverage),
		fmt.Sprintf("%.2f", r.PcIdentity),
		r.Database,
		r.Accession,
		r.Product,
		r.Resistance,
	}
	return strings.Join(fields, sep)
}

// FormatHeader returns the header line.
func FormatHeader(sep string) string {
	return strings.Join(Columns, sep)
}

// SortResults sorts results by sequence name and start position.
func SortResults(results []Result) {
	slices.SortFunc(results, func(a, b Result) int {
		if c := cmp.Compare(a.Sequence, b.Sequence); c != 0 {
			return c
		}
		return cmp.Compare(a.Start, b.Start)
	})
}

// DedupHits removes duplicate hits (same qseqid, qstart, qend).
func DedupHits(hits []blast.Hit) []blast.Hit {
	seen := make(map[string]bool)
	var result []blast.Hit
	for _, h := range hits {
		key := fmt.Sprintf("%s~%d~%d", h.QSeqID, h.QStart, h.QEnd)
		if seen[key] {
			continue
		}
		seen[key] = true
		result = append(result, h)
	}
	return result
}

// ProcessHits processes BLAST hits into Result rows.
func ProcessHits(hits []blast.Hit, fileName, dbName string, minCov float64, noPath bool) []Result {
	hits = DedupHits(hits)
	var results []Result
	for _, hit := range hits {
		r, ok := FormatResult(hit, fileName, dbName, minCov, noPath)
		if ok {
			results = append(results, *r)
		}
	}
	SortResults(results)
	return results
}
