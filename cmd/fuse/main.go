// Fuse: Mass screening of contigs for antimicrobial resistance
// or virulence genes. Go reimplementation of ABRicate.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"strings"
	"sync"

	"github.com/fmj/fuse/internal/blast"
	"github.com/fmj/fuse/internal/database"
	"github.com/fmj/fuse/internal/fasta"
	"github.com/fmj/fuse/internal/report"
	"github.com/fmj/fuse/internal/seqio"
	"github.com/fmj/fuse/internal/summary"
)

// Version of Fuse
const Version = "1.0.0"

func main() {
	// Options
	var (
		debug    bool
		quiet    bool
		check    bool
		setupdb  bool
		list     bool
		summary_ bool
		identity bool
		noheader bool
		csv      bool
		nopath   bool
		version  bool
		threads  int
		jobs     int
		minid    float64
		mincov   float64
		datadir  string
		db       string
		fofn     string
		outfmt   string
	)

	flag.BoolVar(&debug, "debug", false, "Verbose debug output")
	flag.BoolVar(&quiet, "quiet", false, "Quiet mode, no stderr output")
	flag.BoolVar(&check, "check", false, "Check dependencies are installed")
	flag.BoolVar(&setupdb, "setupdb", false, "Format all the BLAST databases")
	flag.BoolVar(&list, "list", false, "List included databases")
	flag.BoolVar(&summary_, "summary", false, "Summarize multiple reports into a table")
	flag.BoolVar(&identity, "identity", false, "Use %IDENTITY instead of %COVERAGE in summary table")
	flag.BoolVar(&noheader, "noheader", false, "Suppress column header row")
	flag.BoolVar(&csv, "csv", false, "Output CSV instead of TSV")
	flag.BoolVar(&nopath, "nopath", false, "Strip filename paths from FILE column")
	flag.BoolVar(&version, "version", false, "Print version and exit")
	flag.IntVar(&threads, "threads", 1, "Use this many BLAST+ threads per file")
	flag.IntVar(&jobs, "jobs", 1, "Run this many files in parallel (increasing --jobs usually speeds up more than --threads)")
	flag.Float64Var(&minid, "minid", 80, "Minimum DNA %identity")
	flag.Float64Var(&mincov, "mincov", 80, "Minimum DNA %coverage")
	flag.StringVar(&datadir, "datadir", defaultDataDir(), "Databases folder")
	flag.StringVar(&db, "db", "ncbi", "Database to use")
	flag.StringVar(&fofn, "fofn", "", "Run on files listed in this file")
	flag.StringVar(&outfmt, "outfmt", "tsv", "Output format: tsv csv bed gff json")
	flag.Usage = printUsage
	flag.Parse()

	if version {
		fmt.Printf("fuse %s\n", Version)
		os.Exit(0)
	}

	// Validate options
	if threads < 1 {
		fmt.Fprintln(os.Stderr, "ERROR: Invalid --threads option, must be 1 or higher")
		os.Exit(1)
	}
	if jobs < 1 {
		fmt.Fprintln(os.Stderr, "ERROR: Invalid --jobs option, must be 1 or higher")
		os.Exit(1)
	}
	if minid <= 0 || minid > 100 {
		fmt.Fprintln(os.Stderr, "ERROR: --minid must be between 1 and 100")
		os.Exit(1)
	}
	if mincov < 0 || mincov > 100 {
		fmt.Fprintln(os.Stderr, "ERROR: --mincov must be between 0 and 100")
		os.Exit(1)
	}

	outfmts := []string{"tsv", "csv", "bed", "gff", "json"}
	if !slices.Contains(outfmts, outfmt) {
		fmt.Fprintf(os.Stderr, "ERROR: Invalid --outfmt '%s' (need: tsv csv bed gff json)\n", outfmt)
		os.Exit(1)
	}

	// Handle --check
	if check {
		msg(quiet, "Checking that BLAST+ tools are installed:")
		requireExe(quiet, "blastn", "makeblastdb", "blastdbcmd", "blastx")
		requireExe(quiet, "gzip", "unzip")
		msg(quiet, "OK.")
		os.Exit(0)
	}

	// Handle --list and --setupdb
	if list || setupdb {
		if setupdb {
			msg(quiet, "Setting up databases in", datadir)
			if err := database.SetupAllDatabases(datadir); err != nil {
				fmt.Fprintf(os.Stderr, "ERROR: %v\n", err)
				os.Exit(1)
			}
		}
		dbs, err := database.ListDatabases(datadir)
		if err != nil {
			fmt.Fprintf(os.Stderr, "ERROR: %v\n", err)
			os.Exit(1)
		}
		sep := "\t"
		if csv {
			sep = ","
		}
		if len(dbs) > 0 {
			fmt.Println(strings.Join([]string{"DATABASE", "SEQUENCES", "DBTYPE", "DATE"}, sep))
			for _, d := range dbs {
				fmt.Println(strings.Join([]string{d.Name, fmt.Sprintf("%d", d.Sequences), d.DBType, d.Date}, sep))
			}
		} else {
			msg(quiet, "No databases found in", datadir)
		}
		os.Exit(0)
	}

	// Handle --summary
	if summary_ {
		args := flag.Args()
		if len(args) < 1 {
			fmt.Fprintf(os.Stderr, "ERROR: --summary needs >=1 fuse output files")
			os.Exit(1)
		}
		if identity {
			msg(quiet, "Using %IDENTITY for the summary table instead of %COVERAGE")
		}
		w := bufio.NewWriter(os.Stdout)
		err := summary.Generate(args, summary.Options{
			UseIdentity: identity,
			NoPath:      nopath,
			CSV:         csv,
		}, w)
		if err != nil {
			fmt.Fprintf(os.Stderr, "ERROR: %v\n", err)
			os.Exit(1)
		}
		w.Flush()
		os.Exit(0)
	}

	// Check dependencies (system gzip/unzip; BLAST+ is required and checked per run)
	// any2fasta is optional - we have built-in FASTA/GenBank/EMBL parsers
	requireExe(quiet, "gzip", "unzip")

	// Collect input files
	args := flag.Args()
	if fofn != "" {
		data, err := os.ReadFile(fofn)
		if err != nil {
			fmt.Fprintf(os.Stderr, "ERROR: Could not open --fofn '%s'\n", fofn)
			os.Exit(1)
		}
		var files []string
		for _, line := range strings.Split(string(data), "\n") {
			line = strings.TrimSpace(line)
			if line != "" {
				files = append(files, line)
			}
		}
		if len(files) == 0 {
			fmt.Fprintf(os.Stderr, "ERROR: No files found in --fofn '%s'\n", fofn)
			os.Exit(1)
		}
		msg(quiet, "Using files from --fofn and ignoring any provided on the command line.")
		args = files
	}

	if len(args) == 0 {
		printUsage()
		os.Exit(1)
	}

	// Check database exists
	datadirAbs, _ := filepath.Abs(datadir)
	dbPath := filepath.Join(datadirAbs, db, "sequences")
	dbInfo, err := database.GetDBInfo(dbPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "ERROR: BLAST database not found: %s\n", dbPath)
		os.Exit(1)
	}
	msg(quiet, fmt.Sprintf("Using %s database %s: %d sequences - %s",
		dbInfo.DBType, db, dbInfo.Sequences, dbInfo.Date))

	// Output separator
	sep := "\t"
	if csv {
		sep = ","
	}

	// Print header
	if !noheader {
		fmt.Println(report.FormatHeader(sep))
	}

	// Process files, up to --jobs in parallel.
	// Results are kept per input index so output stays in input order.
	type fileResult struct {
		lines []string
		err   error
	}
	results := make([]fileResult, len(args))

	// Cap jobs at the number of files to avoid spawning idle workers.
	njobs := min(jobs, len(args))

	var wg sync.WaitGroup
	indexCh := make(chan int)
	for w := 0; w < njobs; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for idx := range indexCh {
				file := args[idx]
				msg(quiet, "Processing:", file)

				if _, err := os.Stat(file); os.IsNotExist(err) {
					results[idx].err = fmt.Errorf("'%s' does not exist, or is unreadable", file)
					continue
				}

				// Convert input to FASTA
				records, err := seqio.AnyToFASTA(file)
				if err != nil {
					results[idx].err = fmt.Errorf("Failed to convert '%s': %v", file, err)
					continue
				}

				// Write FASTA to temp file for BLAST
				tmpFile, err := os.CreateTemp("", "fuse-*.fasta")
				if err != nil {
					results[idx].err = err
					continue
				}
				tmpPath := tmpFile.Name()

				w := fasta.NewWriter(tmpFile)
				if err := w.WriteAll(records); err != nil {
					tmpFile.Close()
					os.Remove(tmpPath)
					results[idx].err = err
					continue
				}
				tmpFile.Close()

				// Run BLAST
				hits, err := blast.Run(tmpPath, blast.Options{
					DBPath:    dbPath,
					DBType:    dbInfo.DBType,
					MinID:     minid,
					Threads:   threads,
					EValue:    "1E-20",
					CullLimit: 1,
					MaxTarget: 10000,
				})
				os.Remove(tmpPath)
				if err != nil {
					results[idx].err = fmt.Errorf("BLAST failed on '%s': %v", file, err)
					continue
				}

				var lines []string
				// Process hits
				hits2 := report.ProcessHits(hits, file, db, mincov, nopath)
				msg(quiet, fmt.Sprintf("Found %d genes in %s", len(hits2), file))

				// Format output lines
				lines = make([]string, 0, len(hits2))
				for i := range hits2 {
					lines = append(lines, report.FormatLine(&hits2[i], sep))
				}
				results[idx] = fileResult{lines: lines}

				msg(quiet, "Done.")
			}
		}()
	}

	// Feed work and wait for all workers.
	for i := range args {
		indexCh <- i
	}
	close(indexCh)
	wg.Wait()

	// Print results in input order and report any errors.
	exitCode := 0
	for i := range results {
		if results[i].err != nil {
			fmt.Fprintf(os.Stderr, "ERROR: %v\n", results[i].err)
			exitCode = 1
			continue
		}
		for _, line := range results[i].lines {
			fmt.Println(line)
		}
	}
	if exitCode != 0 {
		os.Exit(exitCode)
	}
}

func defaultDataDir() string {
	// Try relative to executable first
	exe, err := os.Executable()
	if err == nil {
		dir := filepath.Join(filepath.Dir(exe), "..", "db")
		if _, err := os.Stat(dir); err == nil {
			return dir
		}
	}
	// Fall back to working directory
	return "db"
}

func msg(quiet bool, args ...string) {
	if !quiet {
		fmt.Fprintln(os.Stderr, strings.Join(args, " "))
	}
}

func requireExe(quiet bool, names ...string) {
	for _, name := range names {
		_, err := exec.LookPath(name)
		if err != nil {
			fmt.Fprintf(os.Stderr, "ERROR: Could not find '%s'. Please install it and ensure it is in the PATH.\n", name)
			os.Exit(1)
		}
	}
}

func printUsage() {
	fmt.Fprintf(os.Stderr, "FUSE\n  Find and collate amplicons in assembled contigs\n")
	fmt.Fprintf(os.Stderr, "USAGE\n")
	fmt.Fprintf(os.Stderr, "  %% fuse --list\n")
	fmt.Fprintf(os.Stderr, "  %% fuse [options] <contigs.{fasta,gbk,embl}[.gz] ...> > out.tab\n")
	fmt.Fprintf(os.Stderr, "  %% fuse [options] --fofn fileOfFilenames.txt > out.tab\n")
	fmt.Fprintf(os.Stderr, "  %% fuse --summary <out1.tab> <out2.tab> <out3.tab> ... > summary.tab\n")
	fmt.Fprintf(os.Stderr, "\nTIP\n  For many input files, increasing --jobs (run files in parallel) usually\n  speeds things up more than increasing --threads.\n")
	fmt.Fprintf(os.Stderr, "OPTIONS\n")
	flag.CommandLine.SetOutput(os.Stderr)
	flag.PrintDefaults()
}
