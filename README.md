# Fuse

Mass screening of contigs for antimicrobial resistance or virulence genes.

This is a Go reimplementation of [ABRicate](https://github.com/tseemann/abricate) by Torsten Seemann. The name "Fuse" plays on the idea that while BLAST means "explosion", this tool finds the fuse — locating the genes that could spark an outbreak.

## Features

- Detects acquired antimicrobial resistance genes and virulence factors
- Supports 12 built-in databases: NCBI, CARD, ResFinder, ARG-ANNOT, MEGARES, VFDB, EcOH, PlasmidFinder, Ecoli_VF, VICTORS, BacMet2, UPEC/ExPEC VF
- Accepts FASTA, GenBank, EMBL input (gzip/bzip2 compressed)
- Batch file processing with `--fofn`
- Multiple output formats: TSV, CSV
- Summary mode for presence/absence matrix
- Multi-threaded BLAST search
- Fully compatible with original ABRicate output format
- Uses the system-installed NCBI BLAST+ (blastn, blastx, makeblastdb, blastdbcmd)

## Requirements

- Go 1.21+ (for building from source)
- **NCBI BLAST+** installed and on `PATH` (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`). If missing, `fuse` reports an error and exits. Install from the NCBI FTP site or your package manager (e.g. `brew install blast` / `conda install -c bioconda blast`).
- `gzip`, `unzip` (for compressed files)
- `any2fasta` is optional (built-in FASTA/GenBank/EMBL parsers handle most cases)

## Installation

Installing the pre-compiled binary is the easiest option. For this fuse to work, make sure you have
the NCBI BLAST+ tools on your `PATH` (see [Requirements](#requirements)).

### 1. Binary release (recommended)

Download the archive for your OS and architecture from the
[Releases](https://github.com/famiji/fuse/releases) page. For example, on Linux x86-64:

```bash
# Download and unpack
wget https://github.com/famiji/fuse/releases/download/v1.0.0/fuse_1.0.0_Linux_x86_64.tar.gz
tar xzf fuse_1.0.0_Linux_x86_64.tar.gz

# Put the binaries on your PATH
sudo install fuse fuse-getdb /usr/local/bin/
```

### 2. From source with `go install`

```bash
go install github.com/fmj/fuse/cmd/fuse@v1.0.0
go install github.com/fmj/fuse/cmd/fuse-getdb@v1.0.0
```

This puts `fuse` and `fuse-getdb` into `$(go env GOPATH)/bin`.

### 3. From source with `make`

```bash
make build       # compile both binaries into ./bin/
make install     # copy them into ~/go/bin
```

## Setup

```bash
# Format all BLAST databases
fuse --setupdb --datadir db

# List available databases
fuse --list --datadir db
```

## Usage

```bash
# Run on a single file
fuse --datadir db assembly.fa

# Run on multiple files, 4 at a time (recommended for many inputs)
fuse --datadir db --jobs 4 assembly1.fa assembly2.fa assembly3.fa assembly4.fa

# Use a specific database
fuse --datadir db --db card assembly.fa

# Quiet mode (no stderr messages)
fuse --datadir db --quiet assembly.fa

# CSV output
fuse --datadir db --csv assembly.fa

# File of filenames
fuse --datadir db --fofn filelist.txt

# Summary table
fuse --summary sample1.tab sample2.tab

# Download/update databases
fuse-getdb --dbdir db --db ncbi
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--db` | `ncbi` | Database to use |
| `--minid` | `80` | Minimum DNA %identity |
| `--mincov` | `80` | Minimum DNA %coverage |
| `--threads` | `1` | Number of BLAST+ threads per file |
| `--jobs` | `1` | Number of input files processed in parallel. For many inputs, raising this usually speeds up more than `--threads`. |
| `--fofn` | | File of filenames to process |
| `--csv` | | Output CSV instead of TSV |
| `--noheader` | | Suppress column header |
| `--nopath` | | Strip file paths |
| `--quiet` | | Suppress stderr messages |
| `--list` | | List available databases |
| `--setupdb` | | Format all BLAST databases |
| `--summary` | | Summarize reports into presence/absence matrix |
| `--identity` | | Use %IDENTITY in summary instead of %COVERAGE |
| `--datadir` | `db` | Database directory |
| `--check` | | Check dependencies |

## Output Columns

| Column | Description |
|--------|-------------|
| FILE | Source filename |
| SEQUENCE | Contig/sequence name |
| START | Start coordinate |
| END | End coordinate |
| STRAND | + or - |
| GENE | Gene name |
| COVERAGE | Coverage range |
| COVERAGE_MAP | Visual coverage map |
| GAPS | Gap openings/gaps |
| %COVERAGE | Percentage coverage |
| %IDENTITY | Percentage identity |
| DATABASE | Database name |
| ACCESSION | Accession number |
| PRODUCT | Gene product |
| RESISTANCE | Resistance phenotype |

## Project Structure

```
.
├── cmd/
│   ├── fuse/              # Main CLI entry point (AMR/virulence gene screening)
│   └── fuse-getdb/        # Database downloader CLI
├── internal/
│   ├── blast/             # BLAST+ integration (system binaries)
│   │   ├── blast.go       # BLAST search engine
│   │   └── system.go      # System PATH lookup for BLAST+ tools
│   ├── database/          # Database management
│   ├── fasta/             # FASTA I/O
│   ├── getdb/             # Database downloaders
│   ├── report/            # Result formatting
│   ├── seqio/             # Sequence file I/O (compression, format conversion)
│   └── summary/           # Summary table generation
├── db/                    # Pre-built databases
├── test/                  # Test files
└── Makefile
```

## License

[GPL-2.0](LICENSE).

## Author

Go reimplementation based on [ABRicate](https://github.com/tseemann/abricate) by Torsten Seemann.