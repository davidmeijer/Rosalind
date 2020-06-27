#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving Rosalind exercise Implement
                    MotifEnumeration (ba2a).

Usage:              python3 ba2a.py <Rosalind input file>

Arguments:
<Rosalind input file>   .txt input file containing space separated integers
                        k and d on the first line, followed by a collection
                        of DNA strings on consecutive lines.
"""
from sys import argv

def parse_input(filename):
    """Parses integers k and d and the collection DNA strings from input file.

    Arguments:
    k (int):            k-mer length.
    d (int):            number of mismatches.
    dnas (list):        list of DNA strings.
    """
    with open(filename, 'r') as fo:
        k, d = list(map(int, fo.readline().strip().split()))
        dnas = []
        for line in fo:
            dnas.append(line.strip())

    return k, d, dnas

def motif_enumeration(dnas, k, d):
    """Returns (k,d)-motifs for k-mers with d mismatches from all permutations.

    Arguments:
    dnas (list):        list of DNA strings of nucleotides.
    k (int):            k-mer length.
    d (int):            number of maximum allowed mismatches.

    Returns:
    patterns (set):     permutations that have maximum d mismatches with k-mers
                        present in dnas.
    """
    patterns = set()

    


    return patterns

def get_kmers(sequence, kmer_size):
    """Returns list of all k-mers with length kmer_size from sequence.

    Arguments:
    sequence:   str, DNA sequence.
    kmer_size:  int, length of k-mers to create from sequence.

    Returns:
    kmers:      list of strings, all k-mers made from sequence with length
                kmer_size.
    """
    kmers = []

    for start_char in range(0, len(sequence) - kmer_size + 1):
        kmer = sequence[start_char : start_char + kmer_size]
        kmers.append(kmer)

    return kmers

def hamming(a, b):
    """Returns Hamming distance between two strings.

    Arguments:
    a (str):    DNA string a.
    b (str):    DNA string b.

    Returns:
    hamm (int): Hamming distance between DNA strings a and b.

    Input strings have to be the same length.
    """
    hamm = 0
    for elem in list(zip(a, b)):
        if elem[0] != elem[1]:
            hamm += 1

    return hamm

def main():
    # Parse input k, d, dnas from Rosalind input file:
    k, d, dnas = parse_input(argv[1])

    # Find all (k,d)-motifs in all DNA strings:
    all_patterns = motif_enumeration(dnas, k, d)

    # Print results to stdout:
    print(all_patterns)



if __name__ == '__main__':
    main()
