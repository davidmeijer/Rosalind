#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Protein Translation
http://rosalind.info/problems/ptra/
"""
# Imports:
import argparse
from Bio.Seq import translate

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
    
def parse_input(fn):
    """
    Parses Rosalind assignment specific input file.
    
    Args:
        fn (str): file name of input file.
        
    Returns:
        s (str): DNA sequence.
        p (str): protein sequence.
    """
    with open(fn, 'r') as fo:
        s = fo.readline().strip()
        p = fo.readline().strip()
    
    return s, p

# Main code:
def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    s, p = parse_input(args.input)
    
    for i in range(6):
        table = i + 1
        translation = translate(s, to_stop=True, table=table)
        if translation == p:
            print(table)
    
    
if __name__ == '__main__':
    main()
