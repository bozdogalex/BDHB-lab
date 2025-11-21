# Lab 1 - Databases & GitHub

## Summary

This report documents the completion of Lab 1 exercises, which involved running demonstration scripts for NCBI database queries and implementing a FASTA download tool with GC content analysis. All demo scripts executed successfully, and the exercise `ex01_multifasta_gc.py` was implemented and validated with multiple test cases.

## Demonstration Scripts

### demo01_entrez_brca1.py - NCBI Entrez Query and GenBank Retrieval

This script demonstrates programmatic access to NCBI databases using BioPython's Entrez module. The script searched for BRCA1 sequences in Homo sapiens and retrieved the GenBank record.

**Results:**
- Query: `BRCA1[Gene] AND "Homo sapiens"[Organism]`
- Retrieved record: 3070263245 (Homo sapiens isolate AB17 BRCA1 protein gene, partial cds)
- Sequence length: 350 base pairs
- GC content: 0.423 (42.3%)
- First 50 nucleotides: `CCTGATGGGTTGTGTTTGGTTTCTTTCAGCATGATTTTGAAGTCAGAGGA`

The retrieved sequence represents a partial coding sequence for the BRCA1 protein gene, a known tumor suppressor involved in DNA repair mechanisms.

### demo02_seq_ops.py - Basic Sequence Manipulation

This script demonstrates fundamental sequence operations using BioPython's Seq module. Operations were performed on a 40-nucleotide DNA sequence.

**Input sequence:** `ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG`

**Operations and results:**
- Transcription: `AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG`
- Translation: `MAIVMGR*KGAR*` (13 amino acids with 2 stop codons)
- Reverse complement: `CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT`
- GC content: 0.564 (56.4%)
- ATG start codon positions: [0, 12]

The high GC content (56.4%) and presence of multiple ATG codons are consistent with sequences from coding regions.

### demo03_dbsnp.py - dbSNP Database Query

This script queries the NCBI dbSNP database for single nucleotide polymorphisms associated with BRCA1 in Homo sapiens.

**Results:**
- Query returned 5 SNP entries
- All variants mapped to chromosome 17
- Chromosomal positions ranged from 43126982 to 43127349
- Functional classifications: upstream_transcript_variant, 2KB_upstream_variant, intron_variant

**Example entries:**
- SNP_ID: 2552282559 at position 17:43127349
- SNP_ID: 2552281972 at position 17:43127232

All identified variants are located in regulatory or intronic regions, suggesting potential impact on gene expression regulation rather than direct protein sequence alteration.

## Exercise Implementation

### ex01_multifasta_gc.py - FASTA Download and GC Content Analysis

**Objective:** Implement a command-line tool to download FASTA sequences from NCBI using either direct accession numbers or search queries, then calculate GC content for each sequence.

**Implementation approach:**
The `download_fasta()` function was implemented to handle both direct accession retrieval and query-based searches. The function configures BioPython's Entrez module with user credentials, performs database searches when necessary, retrieves FASTA-formatted sequences via efetch, and writes results to specified output files. GC content calculation uses a robust algorithm that handles non-standard nucleotides.

**Test Case 1: Direct Accession Retrieval**
- Command: `--accession NM_000546`
- Database: nuccore
- Result: 1 record downloaded (NM_000546.6)
- GC content: 0.534 (53.4%)
- Sequence identity: TP53 transcript variant 1, mRNA

**Test Case 2: Query-based Retrieval - Genomic Region**
- Command: `--query "TP53[Gene] AND Homo sapiens[Organism]" --retmax 1`
- Database: nuccore
- Result: 1 record downloaded (NG_017013.2)
- GC content: 0.490 (49.0%)
- Sequence identity: TP53 RefSeqGene on chromosome 17, 32,772 bp genomic region

**Test Case 3: Query-based Retrieval - mRNA Transcripts**
- Command: `--query "TP53[Gene] AND Homo sapiens[Organism] AND mRNA[Filter]" --retmax 2`
- Database: nuccore
- Results: 2 records downloaded
  - NM_001126117.2: GC content = 0.514 (51.4%)
  - NM_001126115.2: GC content = 0.519 (51.9%)

**Analysis and Observations:**

GC content variation was observed across different sequence types for the TP53 gene. Transcript sequences (NM_000546.6, NM_001126117.2, NM_001126115.2) exhibited GC contents ranging from 51.4% to 53.4%, which is above the human genome average of approximately 41%. In contrast, the genomic region (NG_017013.2) showed a lower GC content of 49.0%, likely due to the inclusion of intergenic and intronic sequences.

The elevated GC content in TP53 coding sequences is consistent with its classification as a crucial tumor suppressor gene. Higher GC content in coding regions is associated with increased gene expression regulation and structural stability.

**Technical considerations:**
During testing, it was observed that downloading multiple large genomic records can result in extended processing times due to NCBI API response latency. Using specific search filters (e.g., mRNA[Filter]) significantly improves retrieval performance by limiting results to transcript sequences rather than full genomic regions.

**Code validation:**
The implementation successfully handles:
- Entrez email configuration (required by NCBI)
- Optional API key support for higher rate limits
- Both direct accession and query-based retrieval modes
- Proper error handling for empty search results
- Accurate GC content calculation excluding non-ATGC nucleotides
- File I/O with automatic directory creation

## Deliverables

The following files were created as part of this lab assignment:

- `labs/01_intro&databases/submissions/rbals/ex01_multifasta_gc.py` - Completed exercise implementation
- `labs/01_intro&databases/submissions/rbals/notes.md` - This report
- `data/work/rbals/lab01/*.fa` - Downloaded FASTA sequences (local only, excluded from git repository)

All demonstration scripts were executed successfully, and the exercise implementation passed validation with multiple test cases covering both direct accession retrieval and query-based searches.

## Development Methodology

The exercise code was written manually with inline code completion assistance. The notes were initially drafted by hand and subsequently restructured and refined by GitHub Copilot agent for clarity and scientific formatting. All biological observations and technical findings were cross-validated using ChatGPT to seach for relevant resources.
