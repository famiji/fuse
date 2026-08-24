package database

import (
	"os"
	"testing"
)

func TestParseSeqID(t *testing.T) {
	db, gene, acc, res := ParseSeqID("ncbi~~~blaZ~~~NG_055999.1~~~BETA-LACTAM", "default")
	if db != "ncbi" {
		t.Errorf("expected db 'ncbi', got '%s'", db)
	}
	if gene != "blaZ" {
		t.Errorf("expected gene 'blaZ', got '%s'", gene)
	}
	if acc != "NG_055999.1" {
		t.Errorf("expected acc 'NG_055999.1', got '%s'", acc)
	}
	if res != "BETA-LACTAM" {
		t.Errorf("expected res 'BETA-LACTAM', got '%s'", res)
	}
}

func TestParseSeqIDPartial(t *testing.T) {
	// Only 3 parts: DB~~~GENE~~~ACC
	db, gene, acc, res := ParseSeqID("ncbi~~~blaZ~~~NG_055999.1", "default")
	if db != "ncbi" || gene != "blaZ" || acc != "NG_055999.1" || res != "" {
		t.Errorf("got (%s,%s,%s,%s)", db, gene, acc, res)
	}
}

func TestParseSeqIDNoSeparator(t *testing.T) {
	// No separator at all - should use defaultDB
	db, gene, acc, res := ParseSeqID("plain_id", "mydb")
	if db != "mydb" || gene != "plain_id" || acc != "" || res != "" {
		t.Errorf("got (%s,%s,%s,%s)", db, gene, acc, res)
	}
}

func TestFormatSeqID(t *testing.T) {
	id := FormatSeqID("ncbi", "blaZ", "NG_055999.1", "BETA-LACTAM")
	expected := "ncbi~~~blaZ~~~NG_055999.1~~~BETA-LACTAM"
	if id != expected {
		t.Errorf("expected '%s', got '%s'", expected, id)
	}
}

func TestGuessMolType(t *testing.T) {
	// Create a temp FASTA file with nucleotide sequences
	tmpDir := t.TempDir()
	seqPath := tmpDir + "/test.fa"
	content := ">seq1\nATGCATGCATGCATGC\n>seq2\nGGGCCCATTTGGGCCC\n"
	if err := writeFile(seqPath, content); err != nil {
		t.Fatal(err)
	}
	mt, err := guessMolType(seqPath)
	if err != nil {
		t.Fatal(err)
	}
	if mt != "nucl" {
		t.Errorf("expected 'nucl', got '%s'", mt)
	}
}

func TestGuessMolTypeProt(t *testing.T) {
	tmpDir := t.TempDir()
	seqPath := tmpDir + "/test.fa"
	content := ">seq1\nMQPWWWWLLLLPPPP\n>seq2\nMKLFTDESRVPQWLM\n"
	if err := writeFile(seqPath, content); err != nil {
		t.Fatal(err)
	}
	mt, err := guessMolType(seqPath)
	if err != nil {
		t.Fatal(err)
	}
	if mt != "prot" {
		t.Errorf("expected 'prot', got '%s'", mt)
	}
}

func writeFile(path, content string) error {
	return os.WriteFile(path, []byte(content), 0644)
}
