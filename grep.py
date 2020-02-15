#!/usr/bin/env python3
"""Author: David Meijer"""
import argparse
import string
import copy

def define_arguments():
    """Defines possible command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
                        
    return parser
    
class DNA(object):
    """DNA class."""
    def __init__(self, DNAString):
        self.string = DNAString
    
    def __str__(self):
        return self.string
        
    def __eq__(self, other):
        return self.string == other.string
    
    def __hash__(self):
        return hash(self.string)
    
    def __repr__(self):
        return self.string
    
    def reverse(self):
        """Returns reverse of DNA string."""
        self.string = self.string[::-1]
        
        return self
        
    def complement(self):
        """Returns complement of DNA string."""
        table = self.string.maketrans('ACTG', 'TGAC')
        self.string = self.string.translate(table)
        
        return self
        
    def reverseComplement(self):
        """Returns reverse complement of DNA string."""
        self.reverse()
        self.complement()
        return self
        
def getAdjacencyList(S, kmerSize):
    """Returns De Bruijn graph adjacency list for k-mer size."""
    L = set()
    
    strings = copy.deepcopy(S)
    strings += [s.reverseComplement() for s in S]
    
    for dna in strings:
        str1 = dna.string[:kmerSize]
        for i in range(1, len(dna.string) - (kmerSize - 1)):
            str2 = dna.string[i:i + kmerSize]
            L.add((str1, str2))
            str1 = str2
     
    return L

def makeSuperstring(L):
    """Creates superstring from De Bruijn adjacency list."""
    superstring = []
    
    Lgraph = dict(L)
    
    k = list(Lgraph.keys())[0]
    while len(superstring) < len(L) / 2:
        superstring += Lgraph[k][-1]
        k = Lgraph[k]
        
    return ''.join(superstring)
        
def main():
    """Main code."""
    args = define_arguments().parse_args()
    
    with open(args.input, 'r') as fo: 
        lines = [line.strip() for line in fo]
    
    S = []
    for line in lines:
        S.append(DNA(line))
    
    L = getAdjacencyList(S, 2)
    
    superstring = makeSuperstring(L)
    
    print(superstring) # Same as GASM exercise. Changes not yet incorporated.

if __name__ == '__main__':
    main()
