package fasta

import (
	"bytes"
	"strings"
	"testing"
)

func TestReadAll(t *testing.T) {
	input := ">seq1 description one\nATGCAAAA\nTTTTGGGG\n>seq2 desc two\nCCCCCC\n"
	recs, err := ReadAll(strings.NewReader(input))
	if err != nil {
		t.Fatalf("ReadAll failed: %v", err)
	}
	if len(recs) != 2 {
		t.Fatalf("expected 2 records, got %d", len(recs))
	}
	if recs[0].ID != "seq1" {
		t.Errorf("expected ID 'seq1', got '%s'", recs[0].ID)
	}
	if recs[0].Desc != "description one" {
		t.Errorf("expected Desc 'description one', got '%s'", recs[0].Desc)
	}
	if recs[0].Seq != "ATGCAAAATTTTGGGG" {
		t.Errorf("expected Seq 'ATGCAAAATTTTGGGG', got '%s'", recs[0].Seq)
	}
	if recs[1].ID != "seq2" {
		t.Errorf("expected ID 'seq2', got '%s'", recs[1].ID)
	}
	if recs[1].Seq != "CCCCCC" {
		t.Errorf("expected Seq 'CCCCCC', got '%s'", recs[1].Seq)
	}
}

func TestWriteAll(t *testing.T) {
	recs := []Record{
		{ID: "seq1", Desc: "desc", Seq: "ATGCAAAA"},
		{ID: "seq2", Desc: "", Seq: "CCCCCC"},
	}
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if err := w.WriteAll(recs); err != nil {
		t.Fatalf("WriteAll failed: %v", err)
	}
	expected := ">seq1 desc\nATGCAAAA\n>seq2\nCCCCCC\n"
	if buf.String() != expected {
		t.Errorf("expected '%s', got '%s'", expected, buf.String())
	}
}
