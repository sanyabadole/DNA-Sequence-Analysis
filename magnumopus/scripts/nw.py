#!/usr/bin/env python3

import sys

ALIGNMENT_SCORES = { "match": 1, "mismatch": -1, "gap": -1}

D_FLAGS = {"up": 4, "left": 2, "diag": 1}

def initialize_matrix(num_rows: int, num_columns: int) -> list[list[int]]:
    """ intitializing the matrix by making a 0 filled matrix 
    """
    matrix = []
    for _ in range(num_rows):
        row = []
        for _ in range(num_columns):
            row.append(0)
        matrix.append(row)
    return matrix

def fill_matrix(
    seq_a: str,
    seq_b: str,
    match: int,
    mismatch: int,
    gap: int
) -> tuple[list[list[int]], list[list[int]]]:
    
    #making a matrix of appropriate score and direction
    n_rows = len(seq_a)+1
    n_cols = len(seq_b)+1
    score_mat = initialize_matrix(n_rows, n_cols)
    direct_mat = initialize_matrix(n_rows, n_cols)

    #filling first row and column
    for i in range(n_cols):
        score_mat[0][i] = gap*i
        direct_mat[0][i] = D_FLAGS["left"]

    for i in range(n_rows):
        score_mat[i][0] = gap*i
        direct_mat[i][0] = D_FLAGS["up"]
  
    direct_mat[0][0] = 0

    #filling matrices
    for i in range(1, n_rows):
        for j in range(1, n_cols):
            base_a = seq_a[i-1]
            base_b = seq_b[j-1]
            if base_a == base_b:
                match_mismatch = match
            else:
                match_mismatch = mismatch

            up = score_mat[i-1][j] + gap
            left = score_mat[i][j-1] + gap
            diag = score_mat[i-1][j-1] + match_mismatch

            best_score = max(up, left, diag)
            score_mat[i][j] = best_score

            if up == best_score:
                direct_mat[i][j] += D_FLAGS["up"]
            if left == best_score:
                direct_mat[i][j] += D_FLAGS["left"]
            if diag == best_score:
                direct_mat[i][j] += D_FLAGS["diag"]

    return score_mat, direct_mat


def traceback(
    seq_a: str,
    seq_b: str,
    score_mat: list[list[int]],
    direct_mat: list[list[int]]
) -> tuple[tuple[str, str], int]:
    i = len(seq_a)
    j = len(seq_b)

    aln_a = []
    aln_b = []
    score = score_mat[i][j]

    while i != 0 or j != 0:
        direction_bits = direct_mat[i][j]
        # use bitwise operations to check flags
        if D_FLAGS["diag"] & direction_bits:
            i -= 1
            j -= 1
            aln_a.append(seq_a[i])
            aln_b.append(seq_b[j])
            continue
        if D_FLAGS["up"] & direction_bits:
            i -= 1
            aln_a.append(seq_a[i])
            aln_b.append("-")
            continue
        if D_FLAGS["left"] & direction_bits:
            j -= 1
            aln_a.append("-")
            aln_b.append(seq_b[j])
            continue

    aln_a = "".join(aln_a[::-1])
    aln_b = "".join(aln_b[::-1])

    return (aln_a, aln_b), score


def needleman_wunsch(
    seq_a: str,
    seq_b: str,
    match: int,
    mismatch: int,
    gap: int,
) -> tuple[tuple[str, str], int]:
    
    score_mat, direct_mat = fill_matrix(seq_a, seq_b, match, mismatch, gap)
    alns, score = traceback(seq_a, seq_b, score_mat, direct_mat)
    
    return alns, score


if __name__ == '__main__':
    seq_a = sys.argv[1]
    seq_b = sys.argv[2]
    alns, score = needleman_wunsch(seq_a, seq_b, **ALIGNMENT_SCORES)

    print("\n".join(alns)) 