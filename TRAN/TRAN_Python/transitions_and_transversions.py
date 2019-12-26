#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: Transitions and Transversions
"""
import argparse

def define_arguments():
    """Defines possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')

    return parser

def parse_fasta(path_to_input):
    """Parses fasta file into dictionary.

    Args:
        path_to_input (str): path to fasta style input file.

    Returns:
        fasta_dict (dict): dictionary containing parsed fasta file
        as {header : seq, ...}.

    """
    fasta_dict = {}

    with open(path_to_input, 'r') as in_fo:
        for line in in_fo:
            if line.startswith('>'):
                header = line.strip()[1:]
                fasta_dict[header] = []
                continue
            else:
                fasta_dict[header].append(line.strip())

    for header, seq in fasta_dict.items():
        fasta_dict[header] = ''.join(seq)

    return fasta_dict

def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    fasta_dict = parse_fasta(args.input)

    if len(fasta_dict.keys()) == 2:
        seqs = []
        for header in list(fasta_dict.keys()):
            seqs.append([list(nucl) for nucl in fasta_dict[header]])

        if len(seqs[0]) == len(seqs[1]):
            transitions, transversions = 0, 0
            for pair in map(list.__add__, seqs[0], seqs[1]):
                if 'A' in pair and 'G' in pair or 'C' in pair and 'T' in pair:
                    transitions += 1
                else:
                    if pair[0] != pair[1]:
                        transversions += 1

            print(round(transitions / transversions, 11))


if __name__ == '__main__':
    main()