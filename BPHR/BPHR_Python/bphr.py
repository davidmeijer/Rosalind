#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Base Quality Distribution
http://rosalind.info/problems/bphr/
"""
# Imports:
import argparse
from Bio import SeqIO
import statistics as stat

# Classes and functions:
def define_arguments():
    """
    Defines possible command line arguments.
    
    Returns:
        parser (obj): defines command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
                        
    return parser

def parse_assignment_input(fn):
    """
    Parses assignment specific input.
    
    Args:
        fn (str): file name of Rosalind input file.
        
    Returns:
        q (int): quality threshold value.
        fn_fq (str): file name of FASTQ entries.
    """
    with open(fn, 'r') as fn_fo:
        # Get thresholds from input file:
        q = int(fn_fo.readline().strip())
        
        # Create new file with only FASTQ format sequences:
        fn_fq = ''.join(fn.rsplit('.', 1)[:-1]) + '.fq'
        with open(fn_fq, 'w') as fn_fq_fo:
            for line in fn_fo:
                fn_fq_fo.write(line)
    
    return q, fn_fq

# Main code:
def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    q, fn_fq = parse_assignment_input(args.input)
    
    count = 0
    scores = [record.letter_annotations['phred_quality'] \
              for record in SeqIO.parse(fn_fq, 'fastq')]
    for nucl_scores in zip(*scores):
        if stat.mean(nucl_scores) < q:
            count += 1
    print(count)
    """
    # Group all records together with same header:
    record_dict = {}
    for record in SeqIO.parse(fn_fq, 'fastq'):
        header = record.id
        phred = record.letter_annotations['phred_quality']
        if not header in record_dict:
            record_dict[header] = [phred]
        else:
            record_dict[header].append(phred)
            
    # Calculate if mean quality score exceeds threshold:
    count = 0
    for header,scores in record_dict.items():
        for nucl_scores in zip(*scores):
            if stat.mean(nucl_scores) < q:
                count += 1
                
    print(count)
    """
    
if __name__ == '__main__':
    main()
