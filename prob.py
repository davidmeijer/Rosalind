#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
import math

def define_arguments():
    """Define possible command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Specific Rosalind exercise input file '+\
                        'for "Introduction to Random Strings".')

    return parser

def parse_rosalind_input(path_to_input):
    """Parses Rosalind exercise specific input.

    Args:
        path_to_input (str): path to Rosalind exercise input file.

    Returns:
        seq (str): nucleotide sequence.
        numbers (list): array of numbers (floats)."""
    with open(path_to_input, 'r') as in_fo:
        seq = in_fo.readline().strip()
        numbers = [float(number) for number in in_fo.readline().strip().split()]

    return seq, numbers

def get_nucleotide_count(seq):
    """Get GC and AT count from nucleotide sequence.

    Args:
        seq (str): nucleotide sequence.

    Returns:
        GC_count (int): number of G and C nucleotides in seq.
        AT_count (int): number of A and T nucleotides in seq."""
    GC_count, AT_count = 0, 0

    for nucl in seq:
        if nucl == 'G' or nucl == 'C':
            GC_count += 1
        if nucl == 'A' or nucl == 'T':
            AT_count += 1

    return GC_count, AT_count

def get_log10_probability(GC_count, AT_count, GC_perc):
    """Calculates chance of seq with certain counts with given GC perc.

    Args:
        GC_count (int): number of G and C nucleotides in seq.
        AT_count (int): number of A and T nucleotides in seq.
        GC_perc (float): given GC percentage to construct seq with.

    Returns:


    """
    # Chance on G or C:
    pGC = GC_perc / 2
    # Chance on A or T:
    pAT = (1 - GC_perc) / 2

    log10_probability = math.log10(pGC ** GC_count) + math.log10(pAT ** AT_count)

    return round(log10_probability, 3)

def main():
    args = define_arguments().parse_args()
    seq, GC_percentages = parse_rosalind_input(args.input)

    GC_count, AT_count = get_nucleotide_count(seq)

    # Determine chance of random seq when constructed with certain GC
    # percentage:
    log10_probabilities = []
    for GC_percentage in GC_percentages:
        log10_probability = get_log10_probability(GC_count, AT_count,
                                                  GC_percentage)
        log10_probabilities.append(log10_probability)

    print(' '.join([str(prob) for prob in log10_probabilities]))

if __name__ == '__main__':
    main()
