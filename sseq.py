#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse

def define_arguments():
    """Defines possible command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind exercise specific input file '+\
                        'or relative/absolute path to it.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parses Rosalind exericse specific input file lines.

    Args:
        path_to_input (str): path to specific input file for Rosalind
        exercise.

    Returns:
        fasta_dict (dict): dictionary containing parses fasta input
        lines as {header : s, header : t} with s as nucleotide string
        ant t as a subsequence of a nucleotide string."""
    fasta_dict = {}

    with open(path_to_input, 'r') as in_fo:
        for line in in_fo:
            if line.startswith('>'):
                header = line.strip()[1:]
                fasta_dict[header] = []
                continue
            else:
                fasta_dict[header].append(line.strip())

    for header, seqs in fasta_dict.items():
        fasta_dict[header] = ''.join(seqs)

    return fasta_dict

def find_indices(s, t):
    """Find one collection of indices of t in s.

    Function finds first nucleotide match in sequence, than partitions
    the sequence on that nucleotide and looks for the next nucleotide
    in the final part of the partitioned sequence. For correct indexing,
    the length of the first two partitioned parts are stored.

    Args:
        s (str): nucleotide sequence.
        t (str): nucleotide subsequence.

    Returns:
        indices (list): list of lists with subsequence indices as lists
        of integers of indexes were to find the subsequence in the
        sequence."""
    indices = []

    corr = 1
    for nucl in t:
        index = s.find(nucl) + corr
        indices.append(index)
        red, match, s = s.partition(nucl)
        corr += len(red) + len(match)

    return indices

def main():
    args = define_arguments().parse_args()
    fasta_dict = parse_rosalind_input(args.input)

    s = fasta_dict[list(fasta_dict.keys())[0]]
    t = fasta_dict[list(fasta_dict.keys())[1]]

    # Find indices of t in s:
    indices = find_indices(s, t)

    print(' '.join([str(index) for index in indices]))

if __name__ == '__main__':
    main()
