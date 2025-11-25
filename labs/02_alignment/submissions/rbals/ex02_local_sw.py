#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercise: Local alignment (Smith–Waterman)

Purpose:
  - Load two sequences from the file downloaded in Lab 1 (data/work/<handle>/lab01/).
  - Implement the basic steps of the SW algorithm to obtain a local alignment.

Example run:
  python labs/02_alignment/submissions/rbals/ex02_local_sw.py --fasta data/work/rbals/lab01/test_tp53_mrna.fa --i1 0 --i2 1
"""

import argparse
from pathlib import Path

from Bio import SeqIO

# ===================== Matrix initialization =========================================


def init_score_matrix_local(m: int, n: int):
    """
    Initialize the (m+1) x (n+1) matrix with all values = 0.
    """
    return [[0 for _ in range(n + 1)] for _ in range(m + 1)]


def score_cell_local(
    score, i: int, j: int, a: str, b: str, match: int, mismatch: int, gap: int
):
    """
    Compute the score for a single cell (i, j).
    Steps:
      - diagonal = score[i-1][j-1] + (match if a == b else mismatch)
      - up       = score[i-1][j] + gap
      - left     = score[i][j-1] + gap
    Result = max(0, diagonal, up, left).
    """
    # Match or mismatch score from diagonal
    diagonal = score[i - 1][j - 1] + (match if a == b else mismatch)

    # Gap in sequence 2 (coming from above)
    up = score[i - 1][j] + gap

    # Gap in sequence 1 (coming from left)
    left = score[i][j - 1] + gap

    # Local alignment: score cannot go below 0
    return max(0, diagonal, up, left)


def smith_waterman(seq1: str, seq2: str, match=3, mismatch=-3, gap=-2):
    """Simplified implementation of the Smith–Waterman algorithm."""
    m, n = len(seq1), len(seq2)

    # Matrix initialization
    score = init_score_matrix_local(m, n)

    max_score = 0
    max_pos = (0, 0)

    # Fill the matrix cell by cell
    # + store the highest score and its position
    for i in range(1, m + 1):
        ai = seq1[i - 1]
        for j in range(1, n + 1):
            bj = seq2[j - 1]
            score[i][j] = score_cell_local(score, i, j, ai, bj, match, mismatch, gap)

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    # ================== Backtracking ==================
    # start from the cell with the maximum score and trace back
    align1, align2 = "", ""
    i, j = max_pos
    while i > 0 and j > 0 and score[i][j] > 0:
        # recompute neighbor scores to decide the direction
        diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
        up = score[i - 1][j] + gap
        # Note: left path is the else case below
        _ = score[i][j - 1] + gap  # left score (used in else branch logic)

        if score[i][j] == diag:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif score[i][j] == up:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:  # left
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2, max_score


def calculate_p_distance(align1: str, align2: str) -> float:
    """
    Calculate p-distance from aligned sequences.
    p-distance = proportion of differing positions (excluding gap-only columns).
    """
    if len(align1) != len(align2):
        raise ValueError("Aligned sequences must have the same length")

    differences = 0
    valid_positions = 0

    for a, b in zip(align1, align2):
        # Skip positions where both have gaps
        if a == "-" and b == "-":
            continue
        valid_positions += 1
        if a != b:
            differences += 1

    if valid_positions == 0:
        return 0.0

    return differences / valid_positions


def load_two_sequences(fasta_path: Path, i1: int, i2: int):
    """
    Load the FASTA file and select two sequences by index.
    """
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 2:
        raise SystemExit("[error] The file must contain at least 2 sequences.")
    if not (0 <= i1 < len(recs) and 0 <= i2 < len(recs)):
        raise SystemExit(f"[error] Invalid indices (0..{len(recs) - 1}).")
    return str(recs[i1].seq), str(recs[i2].seq), recs[i1].id, recs[i2].id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--fasta",
        required=True,
        help="Path to your FASTA file in data/work/<handle>/lab01/",
    )
    ap.add_argument(
        "--i1", type=int, default=0, help="Index of the first sequence (default 0)"
    )
    ap.add_argument(
        "--i2", type=int, default=1, help="Index of the second sequence (default 1)"
    )
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[error] File not found: {fasta_path}")

    s1, s2, id1, id2 = load_two_sequences(fasta_path, args.i1, args.i2)

    # Limit to reasonable length for demonstration (still handles different lengths)
    max_len = 500
    if len(s1) > max_len or len(s2) > max_len:
        print(
            f"[info] Sequences are long ({len(s1)}, {len(s2)} bp). Using first {max_len} bp for demonstration."
        )
        s1, s2 = s1[:max_len], s2[:max_len]

    print(f"\nAligning sequences (lengths: {len(s1)}, {len(s2)} bp)...")
    a1, a2, sc = smith_waterman(s1, s2)

    # Calculate p-distance from local alignment
    p_dist = calculate_p_distance(a1, a2) if len(a1) > 0 else 1.0

    print("\n=== Local Alignment (Smith–Waterman) ===")
    print(f"{id1}  vs  {id2}")
    print(f"Alignment length: {len(a1)}")
    print(f"Score: {sc}")
    print(
        f"p-distance (in aligned region): {p_dist:.4f} ({p_dist * 100:.2f}% differences)"
    )
    print("\nLocal alignment:")
    print(a1)
    print(a2)


if __name__ == "__main__":
    main()
