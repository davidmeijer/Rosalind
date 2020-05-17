#!/usr/bin/env python3
"""
Author: David Meijer

Script for solving Rosalind exercise Longest Increasing Subsequence
(lgis).

Usage: "python3 lgis.py <input_file.txt>"
"""

# import statements
from sys import argv

# functions
def parse_input(filename):
	"""Parses positive integer and permutation from input file.
	
	filename: str, file containing input integer and permutation.
	
	Input file consists of two lines: on the first line there is a 
	single integer and on the second line there is a permutation, 
	several numbers seperated each by a space.
	"""
	with open(filename, 'r') as fo:
		n = int(fo.readline().strip())
		perm = list(map(int, fo.readline().strip().split()))
		
	return n, perm
	
def compute_lcs(A, B):
    """Computes Longest Common Subsequences between strings A and B.

    Args:
        A (str): string A.
        B (str): string B.

    Returns:
        X (arr): matrix containing LCS paths for strings A and B."""
    lenA, lenB = len(A), len(B)

    # Create matrix of zeros with dimensions lenA (col), lenB (row):
    X = [[0 for i in range(lenA + 1)] for j in range(lenB + 1)]
    # lenA/lenB + 1 is for inserting a zero column and row!

    # Populate X with scores for calculating LCS:
    for i in range(lenA):
        for j in range(lenB):
            if A[i] == B[j]:
                # +1 if two nucleotides are the same (diagonal):
                X[j + 1][i + 1] = X[j][i] + 1
            else:
                # If not the same, get highest value from top/left:
                X[j + 1][i + 1] = max(X[j + 1][i], X[j][i + 1])

    return X

def backtrack_lcs(matrix, A, B, i, j):
    """Reads out all routes (LCSs) from LCS matrix.

    Args:
        matrix (list of lists): containing LCSs.
        A (str): string A.
        B (str): string B.
        i (int): location in A.
        j (int): location in B.

    Returns:
        lcs (str): longest common substring from matrix."""
    # Return LCS when begin of one of the strings is reached:
    if i == 0 or j == 0:
        return ""

    # If equal return char:
    if A[i - 1] ==  B[j - 1]:
        return backtrack_lcs(matrix, A, B, i - 1, j - 1) + A[i - 1]

    # Choose left over top if larger, otherwise (equal) choose top:
    if matrix[j][i - 1] > matrix[j - 1][i]:
        return backtrack_lcs(matrix, A, B, i - 1, j)

    return backtrack_lcs(matrix, A, B, i, j - 1)
	
	
# main
def main():
	
	# step 1: parse positive integer and permutation from input file
	n, perm = parse_input(argv[1])
	perm = "".join(list(map(str, perm)))
	
	# step 2:
	order = "".join(list(map(str, list(range(1, n + 1)))))
	
	for count in [order, order[::-1]]:
		lcs = compute_lcs(perm, count)
		seq = backtrack_lcs(lcs, perm, count, len(perm), len(count))
		print(" ".join(list(map(str, [char for char in seq]))))
	
if __name__ == "__main__":
	main()
