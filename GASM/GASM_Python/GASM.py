#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: Genome Assembly Using Reads
"""
import argparse
import string

def define_arguments():
    """Defines possible command line arguments.
    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help = 'Rosalind input file.')
    
    return parser
    
def revcomp(dna):
    """Returns reverse complement for DNA string.
    
    Args:
        dna (str): DNA string.
        
    Returns:
        rcdna (str): reverse complement of DNA string.
        
    """
    rcdna  = dna.translate(dna.maketrans('ACTG', 'TGAC'))[::-1]
    
    return rcdna
    
def DeBruijn_adjacency_list(dnas, kmer_size):
    """Creates adjacency list for De Bruijn graph of list of kmers.
    
    Args:
        kmers (list): list of kmers.
        
    Returns:
        adjl (set): adjacency list for De Bruijn graph for kmer list.
    
    """
    adjl = set()
    
    dnas += [revcomp(dna) for dna in dnas]
    
    for dna in dnas:
        kmer1 = dna[:kmer_size]
        for i in range(1, len(dna) - (kmer_size - 1)):
            kmer2 = dna[i:i + kmer_size]
            adjl.add((kmer1, kmer2))
            kmer1 = kmer2
            
    return adjl
    
def make_superstring(adjl):
    """Creates superstring from De Bruijn adjacency list.
    
    Args:
        adjl (set): adjacency list for De Bruijn graph for kmer list.
        
    Returns:
        superstring (str): superstring from De Bruijn adjacency list.
        
    """
    superstring = []
    
    adjg = dict(adjl)
    
    k = list(adjg.keys())[0]
    while len(superstring) < len(adjl) / 2:
        superstring += adjg[k][-1]
        k = adjg[k]
        
    return ''.join(superstring)

def main():
    """Main code.
    
    """
    args = define_arguments().parse_args()
    
    with open(args.input, 'r') as fo:
        dnas = [dna.strip() for dna in fo]
    
    
    for i in range(len(dnas[0]), 0, -1):
        try:
            adjl = DeBruijn_adjacency_list(dnas, i)
            
            superstring = make_superstring(adjl)
    
            print(superstring)
            
            print('kmer size {0} cycles'.format(i))
            
            break
            
        except:
            print('kmer size {0} does not cycle'.format(i))
    
    exit()
    with open('out_gasm.txt', 'w') as fo:
        fo.write(superstring)
    
if __name__ == '__main__':
    main()
