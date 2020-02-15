#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: Constructing a De Bruijn Graph
"""
import argparse
import string
from operator import itemgetter

def define_arguments():
    """Defines possible command line arguments.
    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help = 'Rosalind input file.')
    
    return parser

def revcomp(dna):
    """Creates reverse complement of DNA string.
    
    Args:
        dna (str): DNA string.
        
    Returns:
        rcdna (str): Reverse complement of DNA string.
        
    """
    table = dna.maketrans('ACTG', 'TGAC')

    rcdna = dna.translate(table)[::-1]
    
    return rcdna
    
def DeBruijn_adjacency_list(dnas):
    """Creates a De Bruijn graph adjacency list from a DNA string list.
    
    Args:
        dnas (list): list of DNA strings.
    
    Returns:
        DeBruijn (set): set of De Bruijn adjacencies from DNA strings.
    
    """
    DeBruijn = set()
    
    dnas += [revcomp(dna) for dna in dnas]
    
    for dna in dnas:
        kmers = (dna[:-1], dna[1:])
        DeBruijn.add(kmers)
    
    return sorted(DeBruijn)
    
def main():
    """Main code.
    
    """
    args = define_arguments().parse_args()
    
    with open(args.input, 'r') as fo:
        dnas = [line for line in [line.strip() for line in fo]]
        
    adjl = DeBruijn_adjacency_list(dnas)
    
    with open('out_dbru.txt', 'w') as fo:
        for adj in adjl:
            fo.write('({0}, {1})\n'.format(adj[0], adj[1]))
    
if __name__ == '__main__':
    main()
