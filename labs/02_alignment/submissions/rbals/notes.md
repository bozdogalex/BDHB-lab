# Lab 2 - Sequence Alignment

## Summary

This report documents the completion of Lab 2 exercises, involving the implementation of global (Needleman–Wunsch) and local (Smith–Waterman) alignment algorithms and their application to **TP53** mRNA sequences. Global alignment (end-to-end) and local alignment (best matching subsequences) were both implemented with dynamic programming. We validated the implementations on human TP53 transcript variants and additional sequences for contrast. Key algorithmic choices and when to prefer each method are noted throughout.

## Dataset Used

**Primary sequences:** Two TP53 mRNA transcript variants from *Homo sapiens*
- Sequence 1: NM_001126117.2 (2,063 bp) — TP53 transcript variant 2
- Sequence 2: NM_001126115.2 (2,003 bp) — TP53 transcript variant 4

**Additional comparative sequences:**
- NM_000546.6 — Human TP53 transcript variant 1 (2,512 bp)
- NM_011640.4 — Mouse TP53 (*Mus musculus*) (1,728 bp)
- NM_007294.4 — Human BRCA1 (7,088 bp)

All sequences were downloaded from NCBI (see Lab 1).

## Exercise 1: Global Alignment (Needleman–Wunsch)

**Implementation Details (global):**
1. **Matrix initialization:** (m+1)×(n+1) score matrix; first row/column initialized with cumulative gap penalties.
2. **Recurrence:** for cell (i,j):
   `diag = S[i-1,j-1] + match_or_mismatch`,
   `up = S[i-1,j] + gap`,
   `left = S[i,j-1] + gap`,
   `S[i,j] = max(diag, up, left)`.
3. **Backtracking:** start at bottom-right; move diagonally/up/left to reconstruct the optimal end‑to‑end alignment.

**Test Results:**

*Intra-gene comparison (first 500 bp of TP53 isoforms):*
- Alignment length: 500
- Score: 500 (scheme: match=+1, mismatch=−1, gap=−2)
- **p‑distance:** 0.0000 (0.00% differences)
- *Observation:* Perfect identity in 5′ region reflects shared UTR structure among TP53 splice variants.

*Cross-species comparison (Human vs Mouse TP53):*
- Alignment length: 544
- Score: 62
- **p‑distance:** 0.3621 (36.21% differences)
- *Observation:* Substantial divergence expected between mammalian orthologs; global alignment captures full-length evolutionary changes.

*Cross-gene comparison (Human TP53 vs BRCA1):*
- Alignment length: 528
- Score: -92
- **p‑distance:** 0.5341 (53.41% differences)
- *Observation:* Negative score and high divergence indicate unrelated sequences; global alignment inappropriate for genes without homology.

## Exercise 2: Local Alignment (Smith–Waterman)

**Implementation Details (local):**
1. **Matrix initialization:** all zeros (permits alignment to start anywhere).
2. **Recurrence with floor:** `H[i,j] = max(0, diag, up, left)`.
3. **Backtracking:** start at the max cell; stop when a 0 is reached to return the highest‑scoring local match.

**Test Results:**

*Intra-gene comparison (first 500 bp of TP53 isoforms):*
- Alignment length: 500
- Score: 1500 (scheme: match=+3, mismatch=−3, gap=−2)
- **p‑distance:** 0.0000 (0.00% differences)
- *Observation:* Best local region spans entire window when sequences are highly similar.

*Cross-species comparison (Human vs Mouse TP53):*
- Alignment length: 518
- Score: 692
- **p‑distance:** 0.3031 (30.31% differences)
- *Observation:* Local alignment isolates conserved functional domains; lower p-distance than global (36.21%) indicates stronger conservation in specific regions.

*Cross-gene comparison (Human TP53 vs BRCA1):*
- Best local match found despite overall poor similarity
- *Observation:* Local alignment finds short conserved motifs even in unrelated sequences; more biologically meaningful than forced global alignment.

## Global vs. Local Alignment — When to Use Which

- **Global (Needleman–Wunsch):** sequences of comparable length with expected colinear homology across the full span (orthologs; full transcripts).
- **Local (Smith–Waterman & BLAST):** sequences with conserved subregions/domains or different lengths; also for database search, where BLAST uses a heuristic to approximate optimal local alignment statistics [2,3,10].

## Biological Interpretation

The three-tier comparison strategy reveals key biological insights:

**Intra-gene (TP53 splice variants):** Perfect conservation (0% divergence) in tested regions reflects strong functional constraint on tumor suppressor activity and shared 5′ regulatory elements. Alternative splicing generates isoforms with distinct C-termini while preserving core DNA-binding and transactivation domains [1,8].

**Cross-species (Human vs Mouse TP53):** Moderate divergence (~36% global, ~30% local) illustrates evolutionary conservation under purifying selection. The lower p-distance in local alignment highlights conserved functional motifs (DNA-binding domain, oligomerization region) amid background sequence drift. This pattern is characteristic of orthologous genes with conserved function across mammalian lineages.

**Cross-gene (TP53 vs BRCA1):** High divergence (~53%) and negative global alignment scores confirm these are unrelated genes without significant homology. Even local alignment struggles to find meaningful matches beyond random short motifs. This demonstrates the importance of choosing appropriate alignment strategies—global alignment is unsuitable for sequences lacking colinear homology.

## Technical Considerations

- **Complexity:** naïve dynamic programming for pairwise alignment runs in **O(mn)** time and space for sequences of lengths *m* and *n* [8,9].
- **Affine gaps (Gotoh):** implementing gap‑open vs gap‑extend penalties yields **O(mn)** time with three coupled matrices [4].
- **Linear‑space traceback:** **Hirschberg / Myers–Miller** reduces space to **O(min{m,n})** for optimal alignments [5].
- **Banded DP:** when sequences are similar, confining computation to a diagonal band reduces work to **O(NW)** time (N = shorter length; W = band width) [6].
- **p‑distance:** computed as the proportion of differing sites among aligned (non‑gap‑only) columns; appropriate for closely related sequences [7].

## Implementation Validation

- **Correctness:** initialization, scoring, and traceback verified for both algorithms.
- **Parameterization:** explored match/mismatch/gap schemes; higher match rewards in local alignment helped emphasize truly conserved blocks.
- **Robustness:** p‑distance logic excludes gap‑only columns, properly handling length differences.

## Deliverables

- `labs/02_alignment/submissions/rbals/ex01_global_nw.py` — Global alignment implementation
- `labs/02_alignment/submissions/rbals/ex02_local_sw.py` — Local alignment implementation
- `labs/02_alignment/submissions/rbals/notes.md` — This report

---

## Development Methodology

The exercise code was written manually with inline code completion assistance. The notes were initially drafted by hand and subsequently restructured and refined by GitHub Copilot agent for clarity and scientific formatting. All biological observations and technical findings were cross-validated using ChatGPT to search for relevant resources.

## References

[1] **Needleman & Wunsch (1970)**. *A general method applicable to the search for similarities in the amino acid sequence of two proteins*. J. Mol. Biol. 48:443–453. https://bioinfo.uib.es/~cesc/docencia/doct0708/NeedlemanWunsch.pdf
[2] **Smith & Waterman (1981)**. *Identification of common molecular subsequences*. J. Mol. Biol. 147:195–197. https://courses.cs.duke.edu/compsci260/fall20/resources/AlignmentPapers/1981.smith.waterman.pdf
[3] **Altschul et al. (1990)**. *Basic local alignment search tool (BLAST)*. J. Mol. Biol. 215:403–410. https://pubmed.ncbi.nlm.nih.gov/2231712/
[4] **Gotoh (1982)**. *An improved algorithm for matching biological sequences* (affine gaps). J. Mol. Biol. 162:705–708. https://pubmed.ncbi.nlm.nih.gov/7166760/
[5] **Myers & Miller (1988)**. *Optimal alignments in linear space*. Bioinformatics 4:11–17. https://academic.oup.com/bioinformatics/article/4/1/11/205106
[6] **Chao, Pearson & Miller (1992)**. *Aligning two sequences within a specified diagonal band*. CABIOS 8:481–487. https://pubmed.ncbi.nlm.nih.gov/1422882/
[7] **MEGA Manual** (p‑distance definition). https://www.megasoftware.net/mega1_manual/Distance.html
[8] **LibreTexts** — *The Needleman–Wunsch Algorithm* (complexity O(mn)). https://bio.libretexts.org/Bookshelves/Computational_Biology/Book%3A_Computational_Biology_-_Genomes_Networks_and_Evolution_%28Kellis_et_al.%29/02%3A_Sequence_Alignment_and_Dynamic_Programming/2.05%3A_The_Needleman-Wunsch_Algorithm
[9] **Smith–Waterman (overview)** — local alignment with zero‑floor initialization. https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
[10] **NCBI Site Guide: Sequence Analysis** — BLAST finds regions of local similarity. https://www.ncbi.nlm.nih.gov/guide/sequence-analysis/
