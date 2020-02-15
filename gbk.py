#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
from Bio import Entrez

def define_arguments():
    """
    Defines possible command line arguments.
    
    Returns:
        parser (obj): contains user input command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
                        
    return parser

def main():
    """Main code."""
    args = define_arguments().parse_args()
    
    with open(args.input, 'r') as fo:
        genus, date1, date2 = fo.read().strip().split('\n')
        
    Entrez.email = "david-meijer@live.nl"
    
    search_term = ('"{0}"[Organism]' + 
        ' AND "{1}"[Publication Date] :' +
        ' "{2}"[Publication Date]').format(genus, date1, date2)
    
    handle = Entrez.esearch(
                            db="nucleotide", 
                            term=search_term)
                            
    record = Entrez.read(handle)
    print(record["Count"])
    
if __name__ == '__main__':
    main()
