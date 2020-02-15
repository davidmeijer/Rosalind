#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
from Bio import SeqIO

def define_arguments():
    """Defines possible command line arguments.
    
    Returns:
        parser (obj): command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
    
    return parser
    
def parse_assignment_input(fn):
    """Parses assignment specific input.
    
    Args:
        fn (str): file name of Rosalind input file."""
    with open(fn, 'r') as fn_fo:
        # Get threshold from input file:
        threshold = int(fn_fo.readline().strip())
        
        # Create new file with only FASTQ format sequences:
        fn_fq = ''.join(fn.rsplit('.', 1)[:-1]) + '.fq'
        with open(fn_fq, 'w') as fn_fq_fo:
            for line in fn_fo:
                fn_fq_fo.write(line)
    
    return threshold, fn_fq
    
def main():
    args = define_arguments().parse_args()
    threshold, fn_fq = parse_assignment_input(args.input)
    count = 0
    for record in SeqIO.parse(fn_fq, 'fastq'):
        score = record.letter_annotations['phred_quality']
        if sum(score)/len(score) < threshold:
            count += 1
    print(count)
    
if __name__ == '__main__':
    main()
