#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Finding Genes with ORFs
http://rosalind.info/problems/orfr/
"""
# Imports:
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import translate

# Classes and functions:
def define_arguments():
    """
    Defines possible command line arguments.
    
    Returns:
        parser (obj): possible input command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
                        
    return parser

# Main code:
def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    
    with open(args.input, 'r') as fo:
        dna = fo.readline().strip()
    
    seq = Seq(dna, IUPAC.unambiguous_dna)
    comp = seq.complement()
    revcomp = seq.reverse_complement()
    
    # Messe code from here onwards, clean up later:
    pps = []
    for i in range(3):
        for seq in [comp, revcomp]:
            frame = seq[i:]
            pp = translate(frame)
            pps.append(pp._data)
         
    proteins = []
    for pp in pps:
        pp_slices = pp.split('*')
        protein = []
        orf = False
        for pp_slice in pp_slices:
            for aa in pp_slice:
                if aa == 'M':
                    orf = True
                if orf == True:
                    protein.append(aa)
            if len(protein) != 0:
                proteins.append(''.join(protein))
            protein = []
            orf = False
            
    for prot in proteins:
        print(prot)
    
    for i,length in enumerate([len(prot) for prot in proteins]):
        if length == max([len(prot) for prot in proteins]):
            print(proteins[i])
    
if __name__ == '__main__':
    main()
    
