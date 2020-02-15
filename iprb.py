#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse

def define_arguments():
    """Parses input arguments.

    Function parses input arguments and returns a command line input
    object containing all input arguments.

    Args:
        command line arguments.

    Returns:
        parser (object): parsed command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')
    return parser

def parse_rosalind_input(path_to_input):
    """Parses rosalind input for this specfic exercise.

    Returns:
        k (int): homozygous individuals (AA).
        m (int): heterozygous individuals (Aa).
        n (int): homozygous recessive individuals (aa)."""
    with open(path_to_input, 'r') as in_fo:
        input = in_fo.readline().strip().split(' ')
        input = [int(x) for x in input]
        k, m, n = input[0], input[1], input[2]

    return k, m, n

def chance_dominance_alelle_presentation(k, m, n):
    """Calculates chance dominant allele presentation according Mendel.

    Assumption: any two individuals can mate.

    Args:
        k (int): homozygous individuals (AA).
        m (int): heterozygous individuals (Aa).
        n (int): homozygous recessive individuals (aa).

    Returns:
        p (float): chance to display dominant allele. 5 digits."""
    p = 0

    total = k + m + n

    if k >= 1:
        p += (k / total) * (((k - 1) / (total - 1)) * 1)
        p += (k / total) * ((m / (total - 1)) * 1)
        p += (k / total) * ((n / (total - 1)) * 1)

    if m >= 1:
        p += (m / total) * (((m - 1) / (total - 1)) * (3/4))
        p += (m / total) * ((k / (total - 1)) * 1)
        p += (m / total) * ((n / (total - 1)) * (1/2))

    if n >= 1:
        p += (n / total) * ((k / (total - 1)) * 1)
        p += (n / total) * ((m / (total - 1)) * (1/2))

    return round(float(p), 5)

def main():
    """Main code."""
    args = define_arguments().parse_args()
    k, m, n = parse_rosalind_input(args.input)
    p = chance_dominance_alelle_presentation(k, m, n)
    print(p)

if __name__ == '__main__':
    main()
