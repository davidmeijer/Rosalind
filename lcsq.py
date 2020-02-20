#!/usr/bin/env python3
"""Author: David Meijer"""

import sys
import argparse
import random

sys.setrecursionlimit(1500)

def define_arguments():
    """Defines possible command line arguments.

    Returns:
        parser (object): contains user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')

    return parser

def parse_fasta(fn):
    """Parses Rosalind input.

    Args:
        fn (str): file name of Rosalind input file.
        
    Returns:
        fasta_dict (dict): contains sequenses as {header:sequence,...}."""
    fasta_dict = {}
    
    with open(fn) as fo:
        for line in [line.strip() for line in fo]:
            if line.startswith('>'):
                header = line
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line

    return fasta_dict

def random_string(slen, char_set='ACTG'):
    """Generate a random string of fixed length. 
    
    Default: DNA.
    
    Args:
        slen (int): length of random string.
        char_set (str): char to choose from for random string.
        
    Returns:
        rstring (str): random string of length slen."""
    rstring = ''.join(random.choice(char_set) for i in range(slen))
    
    return rstring
    
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

def main():
    """Main code."""
    args = define_arguments().parse_args()
    fasta = parse_fasta(args.input)
    
    # Get the two sequences from the FASTA dictionary:
    seq1 = fasta[[key for key in fasta.keys()][0]]
    seq2 = fasta[[key for key in fasta.keys()][1]]
    
    # Compute LCS matrix for two strings:
    LCS_matrix = compute_lcs(seq1, seq2)
    
    # Get max LCS length from matrix:
    #print(max([max(x) for x in LCS_matrix]))

    # Read out all LCSs from LCS matrix:
    LCS = backtrack_lcs(LCS_matrix, seq1, seq2, len(seq1), len(seq2))
    print(LCS)

if __name__ == '__main__':
    main()
