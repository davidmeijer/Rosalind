#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def define_arguments():
    """Defines possible command line arguments.
    
    Returns:
        parser (obj): defines command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
                        
    return parser
    
def parse_fasta(fn):
    """Parses FASTA file into dictionary as {header:seq,...}.
    
    Args:
        fn (str): file name of file in FASTA format.
        
    Returns:
        fasta_dict (dict): parsed FASTA file as {header:sequence,...}. """
    fasta_dict = {}
    
    with open(fn, 'r') as fo:
        for line in [line.strip() for line in fo]:
            if line.startswith('>'):
                header = line
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line
    
    return fasta_dict

def main():
    args = define_arguments().parse_args()
    fasta_dict = parse_fasta(args.input)
    
    count = 0
    for header, seq in fasta_dict.items():
        seq_obj = Seq(seq, IUPAC.unambiguous_dna)
        compl = seq_obj.reverse_complement()
        if seq == compl:
            count += 1
            
    print(count)
    
if __name__ == '__main__':
    main()
