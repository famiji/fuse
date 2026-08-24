package getdb

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strings"
)

// DownloadResFinder downloads the ResFinder database.
func DownloadResFinder(dir string, force bool) ([]SeqRecord, error) {
	name := "resfinder_db"
	url := "https://bitbucket.org/genomicepidemiology/" + name + ".git"

	repoDir := filepath.Join(dir, name)
	if _, err := os.Stat(repoDir); err == nil && !force {
		fmt.Fprintf(os.Stderr, "Won't overwrite existing %s (use --force)\n", name)
	} else {
		os.RemoveAll(repoDir)
		fmt.Fprintf(os.Stderr, "Cloning %s to %s\n", url, name)
		if err := gitClone(url, repoDir); err != nil {
			return nil, fmt.Errorf("git clone failed: %w", err)
		}
	}

	// Parse phenotypes.txt
	metaFile := filepath.Join(repoDir, "phenotypes.txt")
	metaData, err := os.ReadFile(metaFile)
	if err != nil {
		return nil, fmt.Errorf("can't read phenotypes.txt: %w", err)
	}

	type anno struct {
		ABX []string
	}
	annoMap := make(map[string]*anno)
	for _, line := range strings.Split(string(metaData), "\n") {
		if strings.HasPrefix(line, "#") || line == "" {
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) < 3 {
			continue
		}
		geneRe := regexp.MustCompile(`^(.*?)_\w+$`)
		m := geneRe.FindStringSubmatch(cols[0])
		if m == nil {
			continue
		}
		gene := m[1]
		if annoMap[gene] == nil {
			annoMap[gene] = &anno{}
		}
		for _, abx := range strings.Split(cols[2], ",") {
			abx = strings.TrimSpace(abx)
			lower := strings.ToLower(abx)
			if strings.Contains(lower, "unknown") || strings.Contains(lower, "notes") || strings.Contains(lower, "none") {
				continue
			}
			annoMap[gene].ABX = append(annoMap[gene].ABX, abx)
		}
	}
	fmt.Fprintf(os.Stderr, "get_resfinder: %s %d genes\n", metaFile, len(annoMap))

	// Read all .fsa files
	pattern := filepath.Join(repoDir, "*.fsa")
	matches, _ := filepath.Glob(pattern)
	var seqs []SeqRecord
	for _, fastaFile := range matches {
		// Fix broken fasta files
		runCommand("sed", "-i.bak", `s/\([A-Z]\)>/\1\n>/gi`, fastaFile)

		recs, err := loadFasta(fastaFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "WARNING: can't read %s: %v\n", fastaFile, err)
			continue
		}
		for _, r := range recs {
			// Parse ID: gene_copy_accession
			re := regexp.MustCompile(`^(.*?)_(\d+)_(\S+)$`)
			m := re.FindStringSubmatch(r.ID)
			if m != nil {
				r.ID = m[1] + "_" + m[2]
				r.ACC = m[3]
			}
			if a, ok := annoMap[strings.Split(r.ID, "_")[0]]; ok {
				r.ABX = append(r.ABX, a.ABX...)
			}
			seqs = append(seqs, r)
		}
	}

	return seqs, nil
}

// DownloadPlasmidFinder downloads the PlasmidFinder database.
func DownloadPlasmidFinder(dir string, force bool) ([]SeqRecord, error) {
	name := "plasmidfinder"
	zipFile := name + ".zip"

	if err := download("https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip", filepath.Join(dir, zipFile), force); err != nil {
		return nil, err
	}
	if err := decompressZip(filepath.Join(dir, zipFile), dir); err != nil {
		return nil, err
	}

	// Find .fsa files
	var seqs []SeqRecord
	pattern := filepath.Join(dir, "*.fsa")
	matches, _ := filepath.Glob(pattern)
	for _, f := range matches {
		recs, err := loadFasta(f)
		if err != nil {
			continue
		}
		for _, r := range recs {
			r.DESC = r.ID
			// Parse: name_accession
			re := regexp.MustCompile(`^(.*)_(([A-Z]+|NC_)\d+(\.\d+)?)$`)
			m := re.FindStringSubmatch(r.ID)
			if m != nil {
				id := strings.TrimRight(m[1], "_")
				r.ID = id
				r.ACC = m[2]
			}
			seqs = append(seqs, r)
		}
	}
	return seqs, nil
}

// DownloadMegaRes downloads the MEGARES database.
func DownloadMegaRes(dir string, force bool) ([]SeqRecord, error) {
	zipFile := "megares.zip"
	if err := download("https://www.meglab.org/downloads/megares_v3.00.zip", filepath.Join(dir, zipFile), force); err != nil {
		return nil, err
	}
	if err := decompressZip(filepath.Join(dir, zipFile), dir); err != nil {
		return nil, err
	}

	// Find megares_drugs_*.fasta
	pattern := filepath.Join(dir, "megares_drugs_*.fasta")
	matches, _ := filepath.Glob(pattern)
	if len(matches) == 0 {
		return nil, fmt.Errorf("no megares fasta files found")
	}

	seqs, err := loadFasta(matches[0])
	if err != nil {
		return nil, err
	}

	var okSeqs []SeqRecord
	for _, s := range seqs {
		// >MEG_372|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_MATE_efflux_pumps|ABEM
		parts := strings.Split(s.ID, "|")
		if len(parts) < 5 {
			continue
		}
		note := ""
		if len(parts) >= 6 {
			note = parts[5]
		}
		if note != "" {
			fmt.Fprintf(os.Stderr, "Skipping %s due to: %s\n", parts[0], note)
			continue
		}
		s.ID = parts[4]
		s.ACC = parts[0]
		s.DESC = strings.Join([]string{parts[1], parts[2], parts[3], parts[4]}, ":")
		okSeqs = append(okSeqs, s)
	}
	return okSeqs, nil
}

// DownloadArgAnnot downloads the ARG-ANNOT database.
func DownloadArgAnnot(dir string, force bool) ([]SeqRecord, error) {
	fastaFile := "arg-annot.fa"
	if err := download(
		"https://www.mediterranee-infection.com/wp-content/uploads/2019/09/ARG-ANNOT_NT_V6_July2019.txt",
		filepath.Join(dir, fastaFile), force); err != nil {
		return nil, err
	}

	// Fix syntax errors: remove backslashes
	data, err := os.ReadFile(filepath.Join(dir, fastaFile))
	if err != nil {
		return nil, err
	}
	data = bytes.ReplaceAll(data, []byte("\\"), []byte(""))
	os.WriteFile(filepath.Join(dir, fastaFile), data, 0644)

	seqs, err := loadFasta(filepath.Join(dir, fastaFile))
	if err != nil {
		return nil, err
	}

	// >(AGly)Aac2-Ie:NC_011896:3039059-3039607:549
	for i := range seqs {
		parts := strings.SplitN(seqs[i].ID, ":", 4)
		if len(parts) >= 3 {
			seqs[i].ID = parts[0]
			seqs[i].ACC = parts[1] + ":" + parts[2]
		}
	}
	return seqs, nil
}

// DownloadCARD downloads the CARD database.
func DownloadCARD(dir string, force bool) ([]SeqRecord, error) {
	tarball := "card.tar.bz2"
	if err := download("https://card.mcmaster.ca/latest/data", filepath.Join(dir, tarball), force); err != nil {
		return nil, err
	}
	if err := extractTarBz2(filepath.Join(dir, tarball), dir); err != nil {
		return nil, err
	}

	jsonFile := filepath.Join(dir, "card.json")
	if _, err := os.Stat(jsonFile); err != nil {
		return nil, fmt.Errorf("could not extract card.json from tarball")
	}

	// Parse CARD JSON
	var card map[string]interface{}
	if err := loadJSON(jsonFile, &card); err != nil {
		return nil, fmt.Errorf("failed to parse card.json: %w", err)
	}

	var seqs []SeqRecord
	for _, v := range card {
		g, ok := v.(map[string]interface{})
		if !ok {
			continue
		}
		modelType, _ := g["model_type"].(string)
		if modelType != "protein homolog model" {
			continue
		}

		id, _ := g["model_name"].(string)
		id = strings.ReplaceAll(id, " ", "_")

		modelSeqs, ok := g["model_sequences"].(map[string]interface{})
		if !ok {
			continue
		}
		seqMap, ok := modelSeqs["sequence"].(map[string]interface{})
		if !ok {
			continue
		}
		// Get first key
		var firstKey string
		for k := range seqMap {
			firstKey = k
			break
		}
		if firstKey == "" {
			continue
		}
		dnaEntry, ok := seqMap[firstKey].(map[string]interface{})
		if !ok {
			continue
		}
		dna, ok := dnaEntry["dna_sequence"].(map[string]interface{})
		if !ok {
			continue
		}

		sequence, _ := dna["sequence"].(string)
		accession, _ := dna["accession"].(string)
		strand, _ := dna["strand"].(string)
		fmin, _ := dna["fmin"].(float64)
		fmax, _ := dna["fmax"].(float64)

		var start, stop float64
		if strand == "-" {
			start = fmax
			stop = fmin
		} else {
			start = fmin
			stop = fmax
		}

		acc := fmt.Sprintf("%s:%d-%d", accession, int(start), int(stop))

		desc, _ := g["ARO_description"].(string)
		if desc == "" {
			desc, _ = g["ARO_accession"].(string)
		}

		// Parse ARO categories for drug classes
		var abx []string
		if cats, ok := g["ARO_category"].(map[string]interface{}); ok {
			for _, cv := range cats {
				cat, ok := cv.(map[string]interface{})
				if !ok {
					continue
				}
				className, _ := cat["category_aro_class_name"].(string)
				if className == "Drug Class" {
					abxName, _ := cat["category_aro_name"].(string)
					abxName = strings.ReplaceAll(abxName, " antibiotic", "")
					abxName = strings.ReplaceAll(abxName, " ", "_")
					abx = append(abx, abxName)
				}
			}
		}

		seqs = append(seqs, SeqRecord{
			ID:   id,
			ACC:  acc,
			DESC: desc,
			SEQ:  strings.ToUpper(sequence),
			ABX:  abx,
		})
	}

	return seqs, nil
}

// DownloadNCBI downloads the NCBI AMRFinderPlus database.
func DownloadNCBI(dir string, force bool) ([]SeqRecord, error) {
	afp := "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
	fastaFile := "amr_cds.ffn"
	tsvFile := "amr_cds.tsv"

	if err := download(afp+"/AMR_CDS.fa", filepath.Join(dir, fastaFile), force); err != nil {
		return nil, err
	}
	if err := download(afp+"/ReferenceGeneCatalog.txt", filepath.Join(dir, tsvFile), force); err != nil {
		return nil, err
	}

	// Load TSV - key by refseq_nucleotide_accession (column 10, 0-indexed)
	tsv, err := loadTabular(filepath.Join(dir, tsvFile), 10)
	if err != nil {
		return nil, err
	}
	fmt.Fprintf(os.Stderr, "[%s] Loaded %d records\n", tsvFile, len(tsv))

	// Parse FASTA
	recs, err := loadFasta(filepath.Join(dir, fastaFile))
	if err != nil {
		return nil, err
	}

	var seqs []SeqRecord
	for _, r := range recs {
		// Parse ID: >WP_061158039.1|NG_050200.1|1|1|blaTEM-156|blaTEM|class_A_beta-lactamase_TEM-156
		parts := strings.Split(r.ID, "|")
		if len(parts) < 7 {
			continue
		}
		pi := parts[0]
		acc := parts[1]
		fp := parts[2]
		fn := parts[3]
		gene := parts[4]
		prod := parts[6]

		// Skip fusion genes
		if fp != "1" || fn != "1" {
			continue
		}

		// Ensure accession has version
		if !strings.Contains(acc, ".") {
			acc += ".1"
		}

		t, ok := tsv[acc]
		if !ok {
			continue
		}

		scope, _ := t["scope"]
		typ, _ := t["type"]
		subtype, _ := t["subtype"]
		if scope != "core" || typ != "AMR" || subtype != "AMR" {
			continue
		}

		prod = strings.ReplaceAll(prod, "_", " ")
		if gene == "" {
			return nil, fmt.Errorf("%s: gene is empty", pi)
		}
		if prod == "" {
			return nil, fmt.Errorf("%s: product is empty", pi)
		}

		refseqAcc, _ := t["refseq_nucleotide_accession"]
		subclass, _ := t["subclass"]
		var abx []string
		if subclass != "" {
			abx = strings.Split(subclass, "/")
		}

		seqs = append(seqs, SeqRecord{
			ID:   gene,
			ACC:  refseqAcc,
			DESC: prod,
			SEQ:  r.SEQ,
			ABX:  abx,
		})
	}

	return seqs, nil
}

// DownloadVFDB downloads the VFDB virulence factor database.
func DownloadVFDB(dir string, force bool) ([]SeqRecord, error) {
	gzFile := "vfdb.fa.gz"
	if err := download("http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz", filepath.Join(dir, gzFile), force); err != nil {
		return nil, err
	}
	if err := decompressGzip(filepath.Join(dir, gzFile), filepath.Join(dir, "vfdb.fa")); err != nil {
		return nil, err
	}

	seqs, err := loadFasta(filepath.Join(dir, "vfdb.fa"))
	if err != nil {
		return nil, err
	}

	// >VFG000676(gb|AAD32411) (lef) anthrax toxin lethal factor precursor [...]
	for i := range seqs {
		re := regexp.MustCompile(`^(\w+)\(\w+\|(\w+)(\.\d+)?\)$`)
		m := re.FindStringSubmatch(seqs[i].ID)
		if m != nil {
			seqs[i].ACC = m[2]
		}
		// Parse description for gene name in parentheses
		descRe := regexp.MustCompile(`^\((.*?)\)`)
		dm := descRe.FindStringSubmatch(seqs[i].DESC)
		if dm != nil {
			seqs[i].ID = dm[1]
		}
	}
	return seqs, nil
}

// DownloadEcoliVF downloads the E. coli virulence factor database.
func DownloadEcoliVF(dir string, force bool) ([]SeqRecord, error) {
	fastaFile := "ecoli_vf.ffn"
	if err := download("https://github.com/phac-nml/ecoli_vf/raw/master/data/repaired_ecoli_vfs_shortnames.ffn",
		filepath.Join(dir, fastaFile), force); err != nil {
		return nil, err
	}

	seqs, err := loadFasta(filepath.Join(dir, fastaFile))
	if err != nil {
		return nil, err
	}

	// >VFG000748(gi:2865308) (espF) EspF [EspF (VF0182)] [Escherichia coli O127:H6 str. E2348/69]
	// >SPG000142 (cvac) Escherichia coli cvi cvaC operon. [X57525 434-745]
	for i := range seqs {
		// Parse ID
		idRe := regexp.MustCompile(`^(\w+)(?:\((.*?)\))?$`)
		m := idRe.FindStringSubmatch(seqs[i].ID)
		if m != nil {
			seqs[i].ID = m[1]
			if m[2] != "" {
				seqs[i].ACC = m[2]
			} else {
				seqs[i].ACC = m[1]
			}
		}
		// Parse description: remove strain name at end
		desc := seqs[i].DESC
		// Remove [strain] at end
		desc = regexp.MustCompile(`\s\[.*?\]$`).ReplaceAllString(desc, "")
		// Parse (gene) name
		descRe := regexp.MustCompile(`^(?:\((.*?)\)\s+)?(.*)$`)
		dm := descRe.FindStringSubmatch(desc)
		if dm != nil {
			if dm[1] != "" {
				seqs[i].ID = dm[1]
			}
			seqs[i].DESC = dm[2]
		}
	}
	return seqs, nil
}

// DownloadUpecExpecVF downloads the UPEC/ExPEC VF database.
func DownloadUpecExpecVF(dir string, force bool) ([]SeqRecord, error) {
	tsvFile := "upec_expec_vf.tsv"
	if err := download("https://raw.githubusercontent.com/FordeGenomics/ST167_Code/refs/heads/main/UPEC-ExPEC_VF/UPEC_ExPEC_VF.tsv",
		filepath.Join(dir, tsvFile), force); err != nil {
		return nil, err
	}

	tsv, err := loadTabular(filepath.Join(dir, tsvFile), 0)
	if err != nil {
		return nil, err
	}

	// Try to find the "Gene name" column
	var seqs []SeqRecord
	for _, row := range tsv {
		geneName, _ := row["Gene name"]
		acc, _ := row["Accession/Source"]
		begin, _ := row["Begin"]
		end, _ := row["End"]
		desc, _ := row["Description"]
		class, _ := row["Class"]
		seq, _ := row["Sequence"]

		if geneName == "" {
			continue
		}

		var abx []string
		if class != "" {
			abx = []string{class}
		}

		seqs = append(seqs, SeqRecord{
			ID:   geneName,
			ACC:  acc + ":" + begin + "-" + end,
			DESC: desc,
			SEQ:  strings.ToUpper(seq),
			ABX:  abx,
		})
	}
	return seqs, nil
}

// DownloadEcOH downloads the EcOH serotyping database.
func DownloadEcOH(dir string, force bool) ([]SeqRecord, error) {
	fastaFile := "EcOH.fa"
	if err := download("https://raw.githubusercontent.com/katholt/srst2/master/data/EcOH.fasta",
		filepath.Join(dir, fastaFile), force); err != nil {
		return nil, err
	}

	seqs, err := loadFasta(filepath.Join(dir, fastaFile))
	if err != nil {
		return nil, err
	}

	// >1__fliC__fliC-H1__1 AB028471.1;flagellin;H1
	for i := range seqs {
		idParts := strings.Split(seqs[i].ID, "__")
		if len(idParts) >= 3 {
			seqs[i].ID = idParts[2]
		}
		descParts := strings.Split(seqs[i].DESC, ";")
		if len(descParts) > 0 {
			seqs[i].ACC = descParts[0]
			if len(descParts) > 1 {
				seqs[i].DESC = strings.Join(descParts[1:], " ")
			}
		}
	}
	return seqs, nil
}

// DownloadBacMet2 downloads the BacMet2 database.
func DownloadBacMet2(dir string, force bool) ([]SeqRecord, error) {
	fastaFile := "bacmet2.fa"
	if err := download("http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta",
		filepath.Join(dir, fastaFile), force); err != nil {
		return nil, err
	}

	seqs, err := loadFasta(filepath.Join(dir, fastaFile))
	if err != nil {
		return nil, err
	}

	// >BAC0098|ctpC|sp|P0A502|CTPC_MYCTU Probable manganese/zinc-exporting
	for i := range seqs {
		parts := strings.Split(seqs[i].ID, "|")
		if len(parts) >= 4 {
			seqs[i].ID = parts[1] + "-" + parts[0]
			seqs[i].ACC = parts[2] + ":" + parts[3]
		}
	}
	return seqs, nil
}

// DownloadVictors downloads the VICTORS virulence factor database.
func DownloadVictors(dir string, force bool) ([]SeqRecord, error) {
	ffnFile := "victors.ffn"
	faaFile := "victors.faa"

	if err := download("http://phidias.us/victors/downloads/gen_downloads.php",
		filepath.Join(dir, ffnFile), force); err != nil {
		return nil, err
	}
	if err := download("http://phidias.us/victors/downloads/gen_downloads_protein.php",
		filepath.Join(dir, faaFile), force); err != nil {
		return nil, err
	}

	// Parse FAA for GI -> ACC, DESC mapping
	giMap := make(map[string]struct {
		ACC  string
		DESC string
	})
	faaData, err := os.ReadFile(filepath.Join(dir, faaFile))
	if err != nil {
		return nil, err
	}
	for _, line := range strings.Split(string(faaData), "\n") {
		if !strings.HasPrefix(line, ">") {
			continue
		}
		// >gi|115534244|ref|YP_783826.1| hypothetical protein pCJ01p4 [Campylobacter jejuni]
		re := regexp.MustCompile(`>gi.(\d+).ref.([^|]+). ([^[]+)`)
		m := re.FindStringSubmatch(line)
		if m != nil {
			giMap[m[1]] = struct {
				ACC  string
				DESC string
			}{m[2], strings.TrimSpace(m[3])}
		}
	}

	seqs, err := loadFasta(filepath.Join(dir, ffnFile))
	if err != nil {
		return nil, err
	}

	// >gi|115534241:2616-3152 Campylobacter jejuni plasmid pCJ01, complete sequence
	for i := range seqs {
		re := regexp.MustCompile(`gi.(\d+):(\d+)-(\d+)`)
		m := re.FindStringSubmatch(seqs[i].ID)
		if m != nil {
			gi := m[1]
			if info, ok := giMap[gi]; ok {
				seqs[i].ACC = info.ACC
				seqs[i].DESC = info.DESC
			} else {
				seqs[i].ACC = fmt.Sprintf("gi|%s:%s-%s", gi, m[2], m[3])
				seqs[i].DESC = "hypothetical protein"
			}
		}
	}
	return seqs, nil
}
