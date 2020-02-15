#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
import itertools
import re

def define_arguments():
    """Defines possible command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')

    return parser


def parse_rosalind_input(path_to_input):
    """Parses specific rosalind input file.

    Args:
        path_to_input (str): relative or absolute path to input file
        containing rosalind input data.

    Returns:
        fasta_dict (dict): dictionary containing input DNA fasta
        sequences."""
    fasta_dict = {}

    with open(path_to_input, 'r') as in_fo:
        for line in in_fo:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                fasta_dict[header] = []
                continue
            else:
                fasta_dict[header].append(line)

    for header, seqs in fasta_dict.items():
        fasta_dict[header] = ''.join(seqs)

    return fasta_dict


def permutate_kmers(seq, kmer_size=4):
    """Permutate kmers of available characters in seq.

    Args:
        seq (str): DNA, RNA or peptide sequence.
        kmer_size (int): size of kmers.

    Returns:
        permutated_kmers (list): all possible permutations of kmers
        from sequence of certain kmer size."""
    kmers = []

    elem_list = []
    for char in seq:
        if not char in elem_list:
            elem_list.append(char)

    # Multiply by kmer_size to include multiple instances of char in
    # same kmer (for example: 'AACT'):
    elem_list *= kmer_size

    elem_list_permutations = []

    # Get permutations:
    for perm_tup in itertools.permutations(elem_list, kmer_size):
        elem_list_permutations.append(''.join(perm_tup))

    # Return sorted, unique permutations:
    return sorted(list(set(elem_list_permutations)))

def count_kmers(seq, kmers):
    """Count every kmer in seq.

    Args:
        seq (str): DNA, RNA or peptide sequence.
        kmers (list): all possible permutations of kmers from sequence
        of certain kmer size.

    Returns:
        kmer_counts (list): list of integers in same order as kmers
        counting the number of instances of kmer in sequence."""
    kmer_count = []

    for kmer in kmers:
        regex_kmer = '(?=({0}))'.format(kmer)
        matches = re.finditer(regex_kmer, seq)
        kmer_count.append(len([match.group(1) for match in matches]))

    return kmer_count

def main():
    """Main code."""
    args = define_arguments().parse_args()
    fasta_dict = parse_rosalind_input(args.input)

    for header, seq in fasta_dict.items():
        kmers = permutate_kmers(seq)
        kmer_count = count_kmers(seq, kmers)

        kmer_count = [str(x) for x in kmer_count]
        print(' '.join(kmer_count))

if __name__ == '__main__':
    main()
