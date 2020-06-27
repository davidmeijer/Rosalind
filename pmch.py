#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving Rosalind exercise Perfect Matchings
                    and RNA Secondary Structures (pmch).

Usage:              python3 pmch.py <input.txt>

Arguments:
    input.txt       Rosalind exercise specific input file in fasta format
                    with only one entry (header and sequence).
"""
from argparse import ArgumentParser

def argument_parser():
    """Parses input arguments."""
    parser = ArgumentParser()
    parser.add_argument(dest = "input", type = str, nargs = "+")
    return parser

def parse_fasta(fi):
    """Generator for entries in fasta file."""
    with open(fi, "r") as fo:
        seq = []
        for line in fo:
            if line.startswith(">"):
                if seq:
                    yield key, "".join(seq)
                    seq = []
                key = line.strip()[1:]
            else:
                seq.append(line.strip())
        yield key, "".join(seq)

def faculty(n):
    """Calculates faculty of positive integer n."""
    if n == 1:
        return n
    return n * faculty(n - 1)

def main():
    """Driver code."""
    args = argument_parser().parse_args()
    _, seq = next(parse_fasta(args.input[0]))
    print(faculty(seq.count('A')) * faculty(seq.count('G')))

if __name__ == "__main__":
    main()
