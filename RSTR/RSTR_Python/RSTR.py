#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: RSTR
http://rosalind.info/problems/rstr/
"""
import argparse
from scipy.special import binom


def define_arguments():
    """
    Define possible command line arguments.

    Returns:
        parser (object): contains command line input arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
    return parser

def parse_input(fn):
    """
    Parses Rosalaind input file.

    Args:
        fn (str): Rosalind input filename.

    Returns:
        N (int): positive integer <= 100000.
        x (float): number x between 0 and 1.
        s (str): DNA string of length at most 10 bp.
    """
    with open(fn) as fo:
        N, x = fo.readline().strip().split()
        s = fo.readline().strip().split()[0]
    return int(N), float(x), s

def check_gc_content(dna):
    """
    Calculates GC ratio of DNA string.

    Args:
        dna (string): DNA sequence.

    Returns:
        GC_ratio (float): GC content of DNA sequence.
    """
    GC_count = dna.count('G') + dna.count('C')
    GC_ratio = round(GC_count / len(dna), 1)
    return GC_ratio

def Bernoulli(N, y, x, s):
    """
    Calculates probability of y successes in x trials.

    Probability of AbAb from AbAb x AbAb is 4/16 or 0.25.
    Probability of non-AbAb is therefore 12/16 or 0.75.

    Args:
        N (int): number of trials.
        y (int): number of successes.
        x (float): GC content.
        s (string): DNA string.

    Returns:
        p (float): chance of y successes in x trials.
    """
    nGC = s.count('G') + s.count('C')
    pGC = x/2

    nAT = s.count('T') + s.count('A')
    pAT = (1-x)/2

    prob = (pAT ** nAT) * (pGC ** nGC)
    p = binom(N, y) * ((prob) ** y) * ((1 - (prob)) ** (N - y))
    return p

def p(N, y, x, s):
    """
    Probability of at least N successes in x trials.

    Args:
        N (int): number of trials.
        y (int): number of successes.
        x (float): GC content.
        s (string): DNA string.

    Returns:
        p_sum (float): chance of at least N successes.
    """
    p_sum = 1 - sum(Bernoulli(N, y, x, s) for y in range(y))
    return p_sum

def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    N, x, s = parse_input(args.input)

    print(round(p(N,1,x,s), 3))

if __name__ == '__main__':
    main()