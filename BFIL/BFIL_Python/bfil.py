#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Base Filtration by Quality
http://rosalind.info/problems/bfil/
"""
# Imports:
import argparse
from Bio import SeqIO
import subprocess

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
    
def run_trimmomatic(q, fn_fq):
    """
    Run trimmomatic.
    
    Args:
        q (int): quality threshold value.
        fn_fq (str): file name of FASTQ entries.
    """
    fn_out = ''.join(fn_fq.rsplit('.', 1)[:-1]) + '.out'
    
    cmd = ('trimmomatic SE -phred33 {0} {1} ' +
           'LEADING:{2} TRAILING:{2}' +
           '').format(fn_fq, fn_out, q)
    subprocess.run(cmd, shell=True, check=True)
    
# Main code:
def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    q, fn_fq = parse_assignment_input(args.input)
    run_trimmomatic(q, fn_fq)

if __name__ == '__main__':
    main()
