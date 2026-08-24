// fuse-getdb: Download and format databases for Fuse.
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"

	"github.com/fmj/fuse/internal/database"
	"github.com/fmj/fuse/internal/getdb"
)

func main() {
	var (
		debug  bool
		outdir string
		db     string
		force  bool
	)

	flag.BoolVar(&debug, "debug", false, "Verbose debug output")
	flag.StringVar(&outdir, "dbdir", defaultDataDir(), "Parent folder")
	flag.StringVar(&db, "db", "", "Database to download")
	flag.BoolVar(&force, "force", false, "Force download even if exists")
	flag.Usage = printUsage
	flag.Parse()

	if db == "" {
		fmt.Fprintf(os.Stderr, "ERROR: Please choose a --db from: %s\n", getdb.AvailableDatabases())
		os.Exit(1)
	}

	if !getdb.IsKnownDB(db) {
		fmt.Fprintf(os.Stderr, "ERROR: Unknown --db '%s', choose from: %s\n", db, getdb.AvailableDatabases())
		os.Exit(1)
	}

	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		fmt.Fprintf(os.Stderr, "ERROR: --dbdir '%s' does not exist\n", outdir)
		os.Exit(1)
	}

	dir := filepath.Join(outdir, db)
	if err := os.MkdirAll(dir, 0755); err != nil {
		fmt.Fprintf(os.Stderr, "ERROR: %v\n", err)
		os.Exit(1)
	}

	fmt.Fprintf(os.Stderr, "Setting up '%s' in '%s'\n", db, dir)

	if err := getdb.DownloadDB(db, dir, force); err != nil {
		fmt.Fprintf(os.Stderr, "ERROR: %v\n", err)
		os.Exit(1)
	}

	// Format BLAST database
	seqPath := filepath.Join(dir, "sequences")
	fmt.Fprintf(os.Stderr, "Formatting BLASTN database: %s\n", seqPath)
	if err := database.MakeBlastDatabase(seqPath, db); err != nil {
		fmt.Fprintf(os.Stderr, "ERROR: %v\n", err)
		os.Exit(1)
	}

	fmt.Fprintf(os.Stderr, "Done.\n")
}

func defaultDataDir() string {
	exe, err := os.Executable()
	if err == nil {
		dir := filepath.Join(filepath.Dir(exe), "..", "db")
		if _, err := os.Stat(dir); err == nil {
			return dir
		}
	}
	return "db"
}

func printUsage() {
	fmt.Fprintf(os.Stderr, "SYNOPSIS\n  Download databases for fuse to use\n")
	fmt.Fprintf(os.Stderr, "USAGE\n  fuse-getdb [options] --db DATABASE\n")
	fmt.Fprintf(os.Stderr, "OPTIONS\n")
	flag.CommandLine.SetOutput(os.Stderr)
	flag.PrintDefaults()
	fmt.Fprintf(os.Stderr, "DATABASES\n  %s\n", getdb.AvailableDatabases())
}
