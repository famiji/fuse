package report

import (
	"testing"

	"github.com/fmj/fuse/internal/blast"
	"github.com/fmj/fuse/internal/database"
)

func TestParseSeqID(t *testing.T) {
	// Test standard format: DB~~~GENE~~~ACC~~~RESISTANCE
	db, gene, acc, res := "ncbi", "blaZ", "NG_055999.1", "BETA-LACTAM"
	id := "ncbi~~~blaZ~~~NG_055999.1~~~BETA-LACTAM"
	gotDB, gotGene, gotAcc, gotRes := database.ParseSeqID(id, "default")
	if gotDB != db || gotGene != gene || gotAcc != acc || gotRes != res {
		t.Errorf("expected (%s,%s,%s,%s), got (%s,%s,%s,%s)",
			db, gene, acc, res, gotDB, gotGene, gotAcc, gotRes)
	}
}

func TestParseSeqIDDefault(t *testing.T) {
	// Test unparseable ID - should use defaultDB
	gotDB, gotGene, gotAcc, gotRes := database.ParseSeqID("simple_id", "mydb")
	if gotDB != "mydb" || gotGene != "simple_id" || gotAcc != "" || gotRes != "" {
		t.Errorf("expected (mydb,simple_id,,), got (%s,%s,%s,%s)",
			gotDB, gotGene, gotAcc, gotRes)
	}
}

func TestMiniMap(t *testing.T) {
	// Full coverage should be all '='
	m := MiniMap(1, 100, 100, 0)
	if m != "===============" {
		t.Errorf("expected 15 '=' chars, got '%s' (len=%d)", m, len(m))
	}
}

func TestMiniMapPartial(t *testing.T) {
	// Partial coverage
	m := MiniMap(1, 50, 100, 0)
	if len(m) != 15 {
		t.Errorf("expected length 15, got %d", len(m))
	}
}

func TestFormatResult(t *testing.T) {
	hit := blast.Hit{
		QSeqID:  "contig1",
		QStart:  100,
		QEnd:    500,
		QLen:    1000,
		SSeqID:  "ncbi~~~blaZ~~~NG_055999.1~~~BETA-LACTAM",
		SStart:  1,
		SEnd:    846,
		SLen:    846,
		SStrand: "plus",
		EValue:  0,
		Length:  846,
		PIdent:  96.81,
		Gaps:    0,
		GapOpen: 0,
		STitle:  "penicillin-hydrolyzing class A beta-lactamase BlaZ",
	}

	r, ok := FormatResult(hit, "test.fa", "ncbi", 80, false)
	if !ok {
		t.Fatal("expected hit to pass filter")
	}
	if r.Gene != "blaZ" {
		t.Errorf("expected Gene 'blaZ', got '%s'", r.Gene)
	}
	if r.Database != "ncbi" {
		t.Errorf("expected Database 'ncbi', got '%s'", r.Database)
	}
	if r.Accession != "NG_055999.1" {
		t.Errorf("expected Accession 'NG_055999.1', got '%s'", r.Accession)
	}
	if r.Resistance != "BETA-LACTAM" {
		t.Errorf("expected Resistance 'BETA-LACTAM', got '%s'", r.Resistance)
	}
	if r.Strand != "+" {
		t.Errorf("expected Strand '+', got '%s'", r.Strand)
	}
	if r.PcCoverage != 100.0 {
		t.Errorf("expected PcCoverage 100.0, got %.2f", r.PcCoverage)
	}
}

func TestFormatResultMinCov(t *testing.T) {
	hit := blast.Hit{
		QSeqID:  "contig1",
		QStart:  100,
		QEnd:    200,
		QLen:    1000,
		SSeqID:  "ncbi~~~blaZ~~~NG_055999.1~~~BETA-LACTAM",
		SStart:  1,
		SEnd:    846,
		SLen:    846,
		SStrand: "plus",
		Length:  100,
		Gaps:    0,
		GapOpen: 0,
		PIdent:  90,
		STitle:  "test",
	}

	// Coverage = 100*100/846 = ~11.8%, should fail with minCov=80
	_, ok := FormatResult(hit, "test.fa", "ncbi", 80, false)
	if ok {
		t.Error("expected hit to fail minCov filter")
	}
}

func TestStrandSymbol(t *testing.T) {
	if strandSymbol("minus") != "-" {
		t.Error("expected '-' for minus strand")
	}
	if strandSymbol("plus") != "+" {
		t.Error("expected '+' for plus strand")
	}
}
