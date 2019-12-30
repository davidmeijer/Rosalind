#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Data Formats
http://rosalind.info/problems/frmt/
"""
# Imports:
import argparse
from Bio import Entrez

# Classes and functions:
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
    
def retrieve_GenBank_record(ids):
    """
    Retrieves GenBank records with GenBank IDs.
    
    Args:
        ids (list): list of GenBank IDs as strings.
        
    Returns:
        records (): retrieved GenBank records.
    """
    Entrez.email  = 'david-meijer@live.nl'
    handle = Entrez.efetch(db='nucleotide', id=ids, rettype='fasta')
    records = handle.read().split('\n\n')

    try:
        for record in records:
            record = record.strip().split('\n')
            header, seq = record[0], ''.join(record[1:])
            # Filter out empty headers and sequences:
            if not header == '' and not seq == '':
                yield header, seq
    except:
        print('Could not retrieve GenBank records!')
        
def shortest_seq_from_dict(fasta_dict):
    """
    Retrieves key for the shortest sequence in a {header:seq,...} dict.
    
    Args:
        fasta_dict (dict): as {header:sequence,...}.
        
    Returns:
        fasta_key (str): key of shortest sequence from fasta_dict.
    """
    shortest_header, shortest_seq = None, float('inf')
    
    for header, seq in fasta_dict.items():
        if len(seq) < shortest_seq:
            shortest_header = header
            shortest_seq = len(seq)
            
    return shortest_header
    
# Main code:
def main():
    """
    Main code.
    """
    # Define command line arguments:
    args = define_arguments().parse_args()
    fn = args.input
    
    # Read IDs from input file:
    with open(fn, 'r') as fo:
        ids = fo.read().strip().split()
    
    # Retrieve records with GenBank IDs:
    fasta_dict = {}
    for header, seq in retrieve_GenBank_record(ids):
        fasta_dict[header] = seq
    
    # Retrieve key for shortest sequence:
    header = shortest_seq_from_dict(fasta_dict)
    
    print(header)
    print(fasta_dict[header])
    
if __name__ == '__main__':
    main()
