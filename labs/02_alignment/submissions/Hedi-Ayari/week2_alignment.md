# Week 2 Sequence Alignment PR Checklist

## General
- [x] I have read the [Contributing Guide](https://github.com/bozdogalex/bioinf-y4-lab/blob/main/CONTRIBUTING.md)
- [x] I have read the [Git Workflow Guide](https://github.com/bozdogalex/bioinf-y4-lab/blob/main/docs/git-workflow.md)
- [x] I understand the [Code of Conduct](https://github.com/bozdogalex/bioinf-y4-lab/blob/main/docs/policies.md)

## Demo Scripts
- [x] I ran demo01_pairwise.py successfully (demonstrates Biopython global/local alignment)
- [x] I ran demo02_distance_matrix.py successfully (computes sequence distances)

## Exercises
- [x] I completed ex02_global_nw.py (Needleman-Wunsch global alignment implementation)
- [x] I completed ex03_local_sw.py (Smith-Waterman local alignment implementation)
- [x] Both exercises use proper matrix initialization and scoring functions
- [x] Both exercises include correct backtracking algorithms

## Submissions
- [x] I created the directory `labs/02_alignment/submissions/<github_handle>/`
- [x] I added `Hedi-Ayari_notes.md` with datasets used and reflection on global vs local alignment
- [x] I added the completed `ex02_global_nw.py` to my submissions directory
- [x] I added the completed `ex03_local_sw.py` to my submissions directory

## Datasets
- [x] I used TP53 sequences from `data/sample/tp53_dna_multi.fasta` for testing
- [x] I understand the difference between global and local alignment algorithms
- [x] I can explain when to use each alignment type

## Pull Request
- [x] My PR title follows the format: `lab2(<handle>): Sequence alignment implementation`
- [x] My PR description includes a summary of implemented algorithms
- [x] I targeted the correct base repository and branch (`bozdogalex/bioinf-y4-lab:main`)
- [x] I have not pushed any data files to the repository (they remain in `data/work/`)

## Learning Outcomes
- [x] I understand global vs local alignment differences
- [x] I can implement basic NW and SW algorithms
- [x] I know when to prefer global over local alignment
- [x] I can interpret alignment results and scores