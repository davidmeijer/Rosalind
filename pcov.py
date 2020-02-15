#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
import numpy as np

def define_arguments():
    """Defines possible command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help = 'Rosalind input file.')
                        
    return parser
    
def make_superstring(adjl):
    """Make superstring from De Bruijn adjacency list.
    
    Args:
        adjl (set): De Bruijn edge adjacency set.
        
    Returns:   
        superstring (str): superstring made from adjacency list."""
    superstring = []
    
    adjg = dict(adjl)
    
    k = list(adjg.keys())[0]
    while len(superstring) < len(adjl):
        superstring += adjg[k][-1]
        k = adjg[k]
    
    return ''.join(superstring)

def main():
    """Main code."""
    args = define_arguments().parse_args()
    
    with open(args.input, 'r') as fo:
        kmers = [kmer.strip() for kmer in fo]
        
    DeBruijn = set()
    for kmer in kmers:
        DeBruijn.add((kmer[:-1], kmer[1:]))
    
    superstring = make_superstring(DeBruijn)
    
    with open('out_pcov.txt', 'w') as fo:
        fo.write(superstring)
    
if __name__ == '__main__':
    main()
