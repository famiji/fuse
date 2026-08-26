// Package getdb provides database download and formatting for Fuse.
package getdb

import (
	"archive/tar"
	"archive/zip"
	"cmp"
	"compress/bzip2"
	"compress/gzip"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"strings"

	"github.com/fmj/fuse/internal/fasta"
)

// SeqRecord represents a database sequence with metadata.
type SeqRecord struct {
	ID   string
	ACC  string
	DESC string
	SEQ  string
	ABX  []string
}

// ABXSep is the separator for resistance phenotypes.
const ABXSep = ";"

// IDSep is the separator in formatted FASTA IDs.
const IDSep = "~~~"

// downloadFunc is a function that downloads and processes a specific database.
type downloadFunc func(dir string, force bool) ([]SeqRecord, error)

// databaseRegistry maps database names to their download functions.
var databaseRegistry = map[string]downloadFunc{
	"resfinder":     DownloadResFinder,
	"plasmidfinder": DownloadPlasmidFinder,
	"megares":       DownloadMegaRes,
	"argannot":      DownloadArgAnnot,
	"card":          DownloadCARD,
	"ncbi":          DownloadNCBI,
	"vfdb":          DownloadVFDB,
	"ecoli_vf":      DownloadEcoliVF,
	"upec_expec_vf": DownloadUpecExpecVF,
	"ecoh":          DownloadEcOH,
	"bacmet2":       DownloadBacMet2,
	"victors":       DownloadVictors,
}

// IsKnownDB returns true if the database name is recognized.
func IsKnownDB(name string) bool {
	_, ok := databaseRegistry[name]
	return ok
}

// AvailableDatabases returns a sorted list of available database names.
func AvailableDatabases() string {
	var names []string
	for name := range databaseRegistry {
		names = append(names, name)
	}
	slices.Sort(names)
	return strings.Join(names, " ")
}

// DownloadDB downloads and formats a database.
func DownloadDB(name, dir string, force bool) error {
	fn, ok := databaseRegistry[name]
	if !ok {
		return fmt.Errorf("unknown database: %s", name)
	}

	srcDir := filepath.Join(dir, "src")
	if err := os.MkdirAll(srcDir, 0755); err != nil {
		return err
	}

	seqs, err := fn(srcDir, force)
	if err != nil {
		return err
	}

	// Deduplicate
	seqs = dedupeSeq(seqs)

	// Sort by ID
	slices.SortFunc(seqs, func(a, b SeqRecord) int {
		return cmp.Compare(a.ID, b.ID)
	})

	// Save FASTA
	seqPath := filepath.Join(dir, "sequences")
	if err := saveFasta(seqPath, seqs, name); err != nil {
		return err
	}

	return nil
}

// download downloads a URL to a file, skipping if already exists and !force.
func download(url, dest string, force bool) error {
	if _, err := os.Stat(dest); err == nil && !force {
		fmt.Fprintf(os.Stderr, "Won't re-download existing %s (use --force)\n", dest)
		return nil
	}

	fmt.Fprintf(os.Stderr, "Downloading: %s\n", url)
	resp, err := http.Get(url)
	if err != nil {
		return fmt.Errorf("failed to download %s: %w", url, err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("HTTP %d for %s", resp.StatusCode, url)
	}

	f, err := os.Create(dest)
	if err != nil {
		return err
	}
	defer f.Close()

	_, err = io.Copy(f, resp.Body)
	if err != nil {
		return err
	}

	info, _ := os.Stat(dest)
	fmt.Fprintf(os.Stderr, "Destination: %s\n", dest)
	fmt.Fprintf(os.Stderr, "Filesize: %d bytes\n", info.Size())
	return nil
}

// dedupeSeq removes duplicate sequences.
func dedupeSeq(seqs []SeqRecord) []SeqRecord {
	seen := make(map[string]string)
	var good []SeqRecord
	for _, s := range seqs {
		if _, exists := seen[s.SEQ]; exists {
			fmt.Fprintf(os.Stderr, "WARNING: duplicate %d bp sequence: %s ~ %s\n",
				len(s.SEQ), s.ID, seen[s.SEQ])
		} else {
			good = append(good, s)
		}
		seen[s.SEQ] = s.ID
	}
	fmt.Fprintf(os.Stderr, "dedupe_seq: read %d / kept %d\n", len(seqs), len(good))
	return good
}

// saveFasta writes sequences in ABRicate format: >DB~~~ID~~~ACC~~~RESISTANCE DESC
func saveFasta(path string, seqs []SeqRecord, dbName string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	w := fasta.NewWriter(f)
	var records []fasta.Record
	for _, s := range seqs {
		abx := ""
		if len(s.ABX) > 0 {
			sorted := slices.Clone(s.ABX)
			slices.Sort(sorted)
			abx = strings.Join(sorted, ABXSep)
			abx = strings.ReplaceAll(abx, " ", "_")
		}
		id := strings.Join([]string{dbName, s.ID, s.ACC, abx}, IDSep)
		desc := s.DESC
		if desc == "" {
			desc = s.ID
		}
		records = append(records, fasta.Record{
			ID:   id,
			Desc: desc,
			Seq:  s.SEQ,
		})
	}
	return w.WriteAll(records)
}

// loadFasta reads a FASTA file and returns SeqRecord slice.
func loadFasta(path string) ([]SeqRecord, error) {
	recs, err := fasta.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var result []SeqRecord
	for _, r := range recs {
		s := strings.ToUpper(r.Seq)
		// Normalize: replace non-AGTC with N for nucleotide
		s = strings.Map(func(c rune) rune {
			if c >= 'A' && c <= 'Z' {
				return c
			}
			return 'N'
		}, s)
		result = append(result, SeqRecord{
			ID:   r.ID,
			DESC: r.Desc,
			SEQ:  s,
		})
	}
	return result, nil
}

// decompressGzip decompresses a .gz file.
func decompressGzip(src, dest string) error {
	f, err := os.Open(src)
	if err != nil {
		return err
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		return err
	}
	defer gz.Close()

	out, err := os.Create(dest)
	if err != nil {
		return err
	}
	defer out.Close()

	_, err = io.Copy(out, gz)
	return err
}

// decompressBzip2 decompresses a .bz2 file.
func decompressBzip2(src, dest string) error {
	f, err := os.Open(src)
	if err != nil {
		return err
	}
	defer f.Close()

	bz := bzip2.NewReader(f)

	out, err := os.Create(dest)
	if err != nil {
		return err
	}
	defer out.Close()

	_, err = io.Copy(out, bz)
	return err
}

// decompressZip decompresses a .zip file to a directory.
func decompressZip(src, destDir string) error {
	r, err := zip.OpenReader(src)
	if err != nil {
		return err
	}
	defer r.Close()

	for _, f := range r.File {
		rc, err := f.Open()
		if err != nil {
			return err
		}

		path := filepath.Join(destDir, f.Name)
		if f.FileInfo().IsDir() {
			os.MkdirAll(path, 0755)
			rc.Close()
			continue
		}

		// Create parent directories
		os.MkdirAll(filepath.Dir(path), 0755)

		out, err := os.Create(path)
		if err != nil {
			rc.Close()
			return err
		}

		_, err = io.Copy(out, rc)
		out.Close()
		rc.Close()
		if err != nil {
			return err
		}
	}
	return nil
}

// extractTarBz2 extracts a .tar.bz2 file.
func extractTarBz2(src, destDir string) error {
	// First decompress bz2
	f, err := os.Open(src)
	if err != nil {
		return err
	}
	defer f.Close()

	bz := bzip2.NewReader(f)
	tr := tar.NewReader(bz)

	for {
		header, err := tr.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		path := filepath.Join(destDir, header.Name)
		switch header.Typeflag {
		case tar.TypeDir:
			os.MkdirAll(path, 0755)
		case tar.TypeReg:
			os.MkdirAll(filepath.Dir(path), 0755)
			out, err := os.Create(path)
			if err != nil {
				return err
			}
			if _, err := io.Copy(out, tr); err != nil {
				out.Close()
				return err
			}
			out.Close()
		}
	}
	return nil
}

// gitClone clones a git repository.
func gitClone(url, dest string) error {
	cmd := exec.Command("git", "clone", "--quiet", url, dest)
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

// runCommand runs a system command.
func runCommand(name string, args ...string) error {
	cmd := exec.Command(name, args...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

// loadTabular reads a TSV file and returns a map keyed by the specified column.
func loadTabular(path string, keyCol int) (map[string]map[string]string, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(string(data), "\n")
	var header []string
	result := make(map[string]map[string]string)
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		cols := strings.Split(line, "\t")
		if header == nil {
			header = cols
			continue
		}
		if len(cols) != len(header) {
			continue
		}
		if keyCol >= len(cols) {
			continue
		}
		key := cols[keyCol]
		if _, exists := result[key]; !exists {
			row := make(map[string]string)
			for i, h := range header {
				if i < len(cols) {
					row[h] = cols[i]
				}
			}
			result[key] = row
		}
	}
	return result, nil
}

// loadJSON reads and parses a JSON file.
func loadJSON(path string, v interface{}) error {
	data, err := os.ReadFile(path)
	if err != nil {
		return err
	}
	return json.Unmarshal(data, v)
}
