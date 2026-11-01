# Lab 01 — Notes (cristian-groza)

Environment: local Docker, image `ghcr.io/bozdogalex/bdhb:base` (Python 3.11 + Biopython).

## Confirmation

I ran all the demo scripts and the exercise:

- [x] `demo01_entrez_brca1.py` — BRCA1 (GenBank) search + download + GC summary
- [x] `demo02_seq_ops.py` — transcription, translation, reverse complement, GC, motif search
- [x] `demo03_dbsnp.py` — dbSNP query for BRCA1 (5 SNP IDs)
- [x] `ex01_multifasta_gc.py` — TODOs completed (Entrez esearch/efetch download, GC on the downloaded file),
  run in both query and accession modes

Downloaded files are saved in `data/work/cristian-groza/lab01/` (git-ignored, not committed).

## Observed results

| Script | Result |
|---|---|
| demo01 | BRCA1 hit ID 3411944935, length 339 bp, **GC fraction = 0.398** |
| demo02 | Protein `MAIVMGR*KGAR*`, **GC fraction = 0.564**, ATG at positions [0, 12] |
| demo03 | 5 SNPs on chr17 (~43.13 Mb), upstream/intron variants |
| ex01 `--accession NM_000546` | `NM_000546.6  GC=0.534` |
| ex01 `--query "TP53[Gene] AND Homo sapiens[Organism]" --retmax 3` | `NG_017013.2  GC=0.490`, `PZ385140.1  GC=0.569`, `PZ385139.1  GC=0.569` |

Observation: the TP53 mRNA (NM_000546, GC = 0.534) is more GC-rich than the TP53 genomic
RefSeqGene region (NG_017013, GC = 0.490), as expected since the genomic region includes the
AT-richer introns.
