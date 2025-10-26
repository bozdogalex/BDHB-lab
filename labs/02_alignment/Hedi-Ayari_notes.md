# Week 2 — Sequence Alignment Notes

## Datasets Used
- Used TP53 sequences from `data/sample/tp53_dna_multi.fasta` for testing the alignment algorithms
- Sequences include human, mouse, and other species TP53 DNA sequences

## Reflection: When is global alignment preferred over local alignment?
Global alignment (Needleman-Wunsch) is preferred when:
- You want to align the entire sequences from end to end
- The sequences are similar in length and you expect them to be homologous across their full length
- You're studying conserved regions across entire genes or proteins
- You need to include all parts of the sequences, even if they have low similarity

Local alignment (Smith-Waterman) is better when:
- You want to find the most similar regions between sequences
- Sequences have different lengths or contain unrelated regions
- You're looking for conserved domains or motifs within larger sequences
- You want to identify the best matching subsequence regardless of position

## Demo Scripts Results
- demo01_pairwise.py: Successfully demonstrated global and local alignments using Biopython on short subsequences
- demo02_distance_matrix.py: Computed Hamming distances and p-distances between all pairs of TP53 sequences

## Exercise Results
- ex01_global_nw.py: Implemented complete Needleman-Wunsch algorithm with proper matrix initialization and scoring
- ex02_local_sw.py: Implemented complete Smith-Waterman algorithm with local alignment logic and max score tracking