// Package database provides antimicrobial resistance and virulence factor
// database management for Fuse.
package database

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"

	"github.com/fmj/fuse/internal/blast"
	"github.com/fmj/fuse/internal/fasta"
)

// IDSeparator is the separator used in database sequence IDs.
const IDSeparator = "~~~"

// DBInfo holds metadata about a formatted database.
type DBInfo struct {
	Name      string
	Sequences int
	DBType    string // "nucl" or "prot"
	Date      string // YYYY-MM-DD
	Path      string // path to the sequences file
}

// SeqRecord represents a database sequence with parsed metadata.
type SeqRecord struct {
	Database    string
	Gene        string
	Accession   string
	Resistance  string // semicolon-separated resistance phenotypes
	Description string
	Sequence    string
}

// ParseSeqID parses a database sequence ID in the format:
//
//	DB~~~GENE~~~ACC~~~RESISTANCE
//
// Returns the database name, gene, accession, and resistance.
// If the ID doesn't contain the separator, sensible defaults are used.
func ParseSeqID(id string, defaultDB string) (db, gene, acc, resistance string) {
	parts := strings.Split(id, IDSeparator)
	if len(parts) >= 4 {
		return parts[0], parts[1], parts[2], parts[3]
	}
	if len(parts) == 3 {
		return parts[0], parts[1], parts[2], ""
	}
	if len(parts) == 2 {
		return parts[0], parts[1], "", ""
	}
	// Can't parse - use the whole thing as gene name
	return defaultDB, id, "", ""
}

// FormatSeqID creates a database sequence ID from components.
func FormatSeqID(db, gene, acc, resistance string) string {
	return strings.Join([]string{db, gene, acc, resistance}, IDSeparator)
}

// LoadDatabase loads a database FASTA file and returns parsed records.
func LoadDatabase(path string) ([]SeqRecord, error) {
	recs, err := fasta.ReadFile(path)
	if err != nil {
		return nil, err
	}
	dbName := filepath.Base(filepath.Dir(path))
	var result []SeqRecord
	for _, r := range recs {
		db, gene, acc, resistance := ParseSeqID(r.ID, dbName)
		result = append(result, SeqRecord{
			Database:    db,
			Gene:        gene,
			Accession:   acc,
			Resistance:  resistance,
			Description: r.Desc,
			Sequence:    r.Seq,
		})
	}
	return result, nil
}

// ListDatabases scans the data directory for available databases.
func ListDatabases(dataDir string) ([]DBInfo, error) {
	entries, err := os.ReadDir(dataDir)
	if err != nil {
		return nil, err
	}
	var dbs []DBInfo
	for _, entry := range entries {
		if !entry.IsDir() {
			continue
		}
		name := entry.Name()
		seqPath := filepath.Join(dataDir, name, "sequences")
		if _, err := os.Stat(seqPath); os.IsNotExist(err) {
			continue
		}
		info, err := GetDBInfo(seqPath)
		if err != nil {
			continue
		}
		info.Name = name
		info.Path = seqPath
		dbs = append(dbs, *info)
	}
	sort.Slice(dbs, func(i, j int) bool {
		return dbs[i].Name < dbs[j].Name
	})
	return dbs, nil
}

// GetDBInfo retrieves database information using blastdbcmd.
func GetDBInfo(seqPath string) (*DBInfo, error) {
	// First check if BLAST index exists
	if !hasBlastIndex(seqPath) {
		return nil, fmt.Errorf("BLAST database not found: %s", seqPath)
	}

	// Use blastdbcmd -info
	blastdbcmd, err := blast.BinaryPath("blastdbcmd")
	if err != nil {
		return nil, err
	}
	cmd := exec.Command(blastdbcmd, "-info", "-db", seqPath)
	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, fmt.Errorf("blastdbcmd failed: %w", err)
	}

	info := &DBInfo{}
	res := string(output)

	// Parse number of sequences
	seqRe := regexp.MustCompile(`([\d,]+)\s+sequences;`)
	m := seqRe.FindStringSubmatch(res)
	if m != nil {
		n, _ := strconv.Atoi(strings.ReplaceAll(m[1], ",", ""))
		info.Sequences = n
	}

	// Parse database type
	if strings.Contains(res, "total residues") {
		info.DBType = "prot"
	} else if strings.Contains(res, "total bases") {
		info.DBType = "nucl"
	}

	// Parse date
	dateRe := regexp.MustCompile(`Date:\s+(\w+)\s+(\d+),\s+(\d{4})`)
	dm := dateRe.FindStringSubmatch(res)
	if dm != nil {
		info.Date = fmt.Sprintf("%s-%s-%s", dm[3], dm[1], dm[2])
	}

	return info, nil
}

// hasBlastIndex checks if BLAST index files exist for the given path.
func hasBlastIndex(seqPath string) bool {
	// Check for nucleotide index (.nhr, .nin, .nsq) or protein index (.phr, .pin, .psq)
	if fileExists(seqPath+".nhr") || fileExists(seqPath+".nin") {
		return true
	}
	if fileExists(seqPath+".phr") || fileExists(seqPath+".pin") {
		return true
	}
	return false
}

func fileExists(path string) bool {
	_, err := os.Stat(path)
	return err == nil
}

// MakeBlastDatabase creates a BLAST database from a FASTA file.
func MakeBlastDatabase(seqPath, title string) error {
	// Determine type (nucl or prot) by checking sequence content
	dbtype, err := guessMolType(seqPath)
	if err != nil {
		return err
	}

	// Remove old index files
	pattern := seqPath + ".[np]??"
	matches, _ := filepath.Glob(pattern)
	for _, m := range matches {
		os.Remove(m)
	}
	// Also remove .n?? files
	matches2, _ := filepath.Glob(seqPath + ".n??")
	for _, m := range matches2 {
		os.Remove(m)
	}

	makeblastdb, err := blast.BinaryPath("makeblastdb")
	if err != nil {
		return err
	}
	cmd := exec.Command(makeblastdb, "-in", seqPath, "-title", title, "-dbtype", dbtype, "-logfile", os.DevNull)
	if output, err := cmd.CombinedOutput(); err != nil {
		return fmt.Errorf("makeblastdb failed: %w: %s", err, string(output))
	}
	return nil
}

// guessMolType determines if a FASTA file contains nucleotide or protein sequences.
func guessMolType(seqPath string) (string, error) {
	recs, err := fasta.ReadFile(seqPath)
	if err != nil {
		return "", err
	}
	if len(recs) == 0 {
		return "nucl", nil
	}

	var letters strings.Builder
	for _, r := range recs {
		letters.WriteString(r.Seq)
	}
	s := letters.String()
	total := len(s)
	if total == 0 {
		return "nucl", nil
	}

	// Count non-AGTC characters
	nonAGTC := 0
	for _, c := range s {
		switch c {
		case 'A', 'G', 'T', 'C', 'a', 'g', 't', 'c':
		default:
			nonAGTC++
		}
	}

	if float64(nonAGTC) > 0.5*float64(total) {
		return "prot", nil
	}
	return "nucl", nil
}

// SetupAllDatabases formats all databases in the data directory.
func SetupAllDatabases(dataDir string) error {
	entries, err := os.ReadDir(dataDir)
	if err != nil {
		return err
	}
	for _, entry := range entries {
		if !entry.IsDir() {
			continue
		}
		name := entry.Name()
		seqPath := filepath.Join(dataDir, name, "sequences")
		if _, err := os.Stat(seqPath); os.IsNotExist(err) {
			continue
		}
		fmt.Fprintf(os.Stderr, "Setting up database: %s\n", name)
		if err := MakeBlastDatabase(seqPath, name); err != nil {
			return fmt.Errorf("failed to setup %s: %w", name, err)
		}
	}
	return nil
}

// DefaultDataDir returns the default database directory relative to the executable.
func DefaultDataDir() string {
	exe, err := os.Executable()
	if err != nil {
		return filepath.Join("..", "db")
	}
	return filepath.Join(filepath.Dir(exe), "..", "db")
}
