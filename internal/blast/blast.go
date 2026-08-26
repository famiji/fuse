// Package blast provides BLAST+ integration for sequence alignment.
package blast

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"strconv"
	"strings"
)

// Hit represents a single BLAST hit with all fields used by Fuse.
type Hit struct {
	QSeqID  string
	QStart  int
	QEnd    int
	QLen    int
	SSeqID  string
	SStart  int
	SEnd    int
	SLen    int
	SStrand string
	EValue  float64
	Length  int
	PIdent  float64
	Gaps    int
	GapOpen int
	STitle  string
}

// BlastFields are the output fields requested from BLAST.
var BlastFields = []string{
	"qseqid", "qstart", "qend", "qlen",
	"sseqid", "sstart", "send", "slen", "sstrand",
	"evalue", "length", "pident", "gaps", "gapopen",
	"stitle",
}

// Options configures the BLAST search.
type Options struct {
	DBPath    string
	DBType    string // "nucl" or "prot"
	MinID     float64
	Threads   int
	EValue    string
	CullLimit int
	MaxTarget int
}

// buildArgs constructs the BLAST command-line arguments. The query is passed
// either via -query inputFile, or (when queryPath is empty) via stdin.
func buildArgs(opts Options, queryPath string) (string, []string) {
	format := "6 " + strings.Join(BlastFields, " ")
	args := []string{
		"-db", opts.DBPath,
		"-outfmt", format,
		"-num_threads", strconv.Itoa(opts.Threads),
		"-evalue", opts.EValue,
		"-culling_limit", strconv.Itoa(opts.CullLimit),
		"-max_target_seqs", strconv.Itoa(opts.MaxTarget),
	}
	if queryPath != "" {
		args = append(args, "-query", queryPath)
	}

	if opts.DBType == "nucl" {
		args = append(args, "-task", "blastn", "-dust", "no",
			"-perc_identity", fmt.Sprintf("%.0f", opts.MinID))
		return "blastn", args
	}
	args = append(args, "-task", "blastx-fast", "-seg", "no")
	return "blastx", args
}

// parseOutput scans BLAST tabular output from r into a slice of Hit.
func parseOutput(r io.Reader) ([]Hit, error) {
	var hits []Hit
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 10*1024*1024)
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) != len(BlastFields) {
			continue
		}

		hit := Hit{}
		hit.QSeqID = fields[0]
		hit.QStart, _ = strconv.Atoi(fields[1])
		hit.QEnd, _ = strconv.Atoi(fields[2])
		hit.QLen, _ = strconv.Atoi(fields[3])
		hit.SSeqID = fields[4]
		hit.SStart, _ = strconv.Atoi(fields[5])
		hit.SEnd, _ = strconv.Atoi(fields[6])
		hit.SLen, _ = strconv.Atoi(fields[7])
		hit.SStrand = fields[8]
		hit.EValue, _ = strconv.ParseFloat(fields[9], 64)
		hit.Length, _ = strconv.Atoi(fields[10])
		hit.PIdent, _ = strconv.ParseFloat(fields[11], 64)
		hit.Gaps, _ = strconv.Atoi(fields[12])
		hit.GapOpen, _ = strconv.Atoi(fields[13])
		hit.STitle = fields[14]

		// Normalize strand: if minus, swap sstart/send
		if hit.SStrand == "minus" {
			hit.SStart, hit.SEnd = hit.SEnd, hit.SStart
		}

		hits = append(hits, hit)
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading BLAST output: %w", err)
	}
	return hits, nil
}

// Run executes BLAST on the input FASTA file (via stdin) and returns parsed hits.
func Run(inputFile string, opts Options) ([]Hit, error) {
	blastCmd, blastArgs := buildArgs(opts, inputFile)
	binPath, err := BinaryPath(blastCmd)
	if err != nil {
		return nil, err
	}
	cmd := exec.Command(binPath, blastArgs...)
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		return nil, fmt.Errorf("failed to create stdout pipe: %w", err)
	}
	cmd.Stderr = os.Stderr

	if err := cmd.Start(); err != nil {
		return nil, fmt.Errorf("failed to start %s: %w", blastCmd, err)
	}

	hits, err := parseOutput(stdout)
	if err != nil {
		cmd.Wait()
		return nil, err
	}

	if err := cmd.Wait(); err != nil {
		return nil, fmt.Errorf("BLAST failed: %w", err)
	}

	return hits, nil
}

// RunWithStdin executes BLAST with input piped via stdin.
func RunWithStdin(fastaData string, opts Options) ([]Hit, error) {
	blastCmd, blastArgs := buildArgs(opts, "")
	binPath, err := BinaryPath(blastCmd)
	if err != nil {
		return nil, err
	}
	cmd := exec.Command(binPath, blastArgs...)
	cmd.Stdin = strings.NewReader(fastaData)
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		return nil, fmt.Errorf("failed to create stdout pipe: %w", err)
	}
	cmd.Stderr = os.Stderr

	if err := cmd.Start(); err != nil {
		return nil, fmt.Errorf("failed to start %s: %w", blastCmd, err)
	}

	hits, err := parseOutput(stdout)
	if err != nil {
		cmd.Wait()
		return nil, err
	}

	if err := cmd.Wait(); err != nil {
		return nil, fmt.Errorf("BLAST failed: %w", err)
	}

	return hits, nil
}
