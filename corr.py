#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
import copy

def define_arguments():
    """Defines possible command line arguments.

    Returns:
        parser (object): contains command line arguments."""
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file filename.')
                        
    return parser

def parse_input(fn):
    """Parses Rosalind input file for this assignment.

    Args:
        fn (str): Rosalind input file filename.

    Returns:
        input (dict): input fasta sequences as {header:seq, ...}."""
    input = {}
    with open(fn) as fo:
        for line in [line.strip() for line in fo]:
            if line.startswith('>'):
                header = line[1:]
                input[header] = []
            else:
                input[header].append(line)
    input.update((header, ''.join(seq)) for header, seq in input.items())
    
    return input

def revcomp(dna):
    """Gets reverse complement of DNA strand.

    Args:
        dna (string): DNA sequence.

    Returns:
        revcomp_dna (string): reverse complement of DNA sequence."""
    rev_dna = dna[::-1]
    td = rev_dna.maketrans('ACTG', 'TGAC')
    revcomp_dna = rev_dna.translate(td)
    
    return revcomp_dna

def main():
    """Main code."""
    args = define_arguments().parse_args()
    input = parse_input(args.input)
    seqs = [seq for header, seq in input.items()]
    # Get loners:
    l = [seq for seq in seqs if seqs.count(seq) == 1]
    # Get reverse complement of loners:
    lrc = [revcomp(seq) for seq in l]
    # Filter loners on reverse complements:
    fl = []
    for i in range(len(l)):
        c_lrc = copy.deepcopy(lrc)
        del c_lrc[i]
        if l[i] not in c_lrc:
            fl.append(l[i])
    print(fl)

if __name__ == '__main__':
    main()
