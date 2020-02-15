#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Read Filtration by Quality
http://rosalind.info/problems/filt/
"""
# Imports:
import argparse
import subprocess
from Bio import SeqIO

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
        p (int): percentage of bases.
        fn_fq (str): file name of FASTQ entries.
    """
    with open(fn, 'r') as fn_fo:
        # Get thresholds from input file:
        q, p = map(int, fn_fo.readline().strip().split())
        
        # Create new file with only FASTQ format sequences:
        fn_fq = ''.join(fn.rsplit('.', 1)[:-1]) + '.fq'
        with open(fn_fq, 'w') as fn_fq_fo:
            for line in fn_fo:
                fn_fq_fo.write(line)
    
    return q, p, fn_fq
    
def fastx_quality_filter(fn, q, p, fn_out):
    """
    Runs FASTX Quality Filter on FASTQ file.
    
    Args:
        fn (str): file name of FASTQ entries.
        q (int): quality threshold value.
        p (int): percentage of bases.
    """
    cmd = ('fastq_quality_filter -q {0} -p {1} -i {2} -o {3}')\
           .format(q, p, fn, fn_out)
    subprocess.run(cmd, shell=True, check=True)
    
    return fn_out
    
# Main code:
def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    q, p, fn_fq = parse_assignment_input(args.input)
    fn_out = ''.join(fn_fq.rsplit('.', 1)[:-1]) + '.out'
    fn_out = fastx_quality_filter(fn_fq, q, p, fn_out)
    print(len([record for record in SeqIO.parse(fn_out, 'fastq')]))
    
if __name__ == '__main__':
    main()
