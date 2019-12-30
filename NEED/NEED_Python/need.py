#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Pairwise Global Alignment
http://rosalind.info/problems/need/
"""
# Imports:
import argparse
from Bio import Entrez
import subprocess

# Classes and functions:
def define_arguments():
    """
    Defines command line arguments.
    
    Returns:
        parser (obj): defines command line arguments.
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
        
def run_needle(seqA, seqB):
    """
    Run alignment with Needle from EMBOSS on two sequences.
    
    Args:
        seqA (str): (DNA) sequence.
        seqB (str): (DNA) sequence.
        
    Returns:
        out (str): file name of output alignment file.
    """
    fn1, fn2, out = './seq1.fasta', './seq2.fasta', '.alignment.txt'
    
    with open(fn1, 'w') as fo1:
        fo1.write('>seq1\n{0}\n'.format(seqA))
    
    with open(fn2, 'w') as fo2:
        fo2.write('>seq2\n{0}\n'.format(seqB))
    
    cmd = ('/Users/david/miniconda3/bin/needle {0} {1}' +
    ' -gapopen {2} -gapextend {3} -outfile {4}' +
    ' -endweight Yes -endopen 10 -endextend 1').format(fn1, fn2, 10, 1, out)

    result = subprocess.run(cmd, shell=True, check=True)

    for fn in [fn1, fn2]:
        subprocess.run('rm {0}'.format(fn), shell=True, check=True)
    
    #print(result.returncode)
    
    return out
    
def parse_needle(fn):
    """
    Parses alignment file and returns alignment score.
    
    Args:
        fn (str): filename of Needle output alignment file.
        
    Returns:
        score (int): alignment score from Needle alignment file.
    """
    score = 0
    
    with open(fn, 'r') as fo:
        for line in [line.strip() for line in fo]:
            if line.startswith('# Score: '):
                score = line.split()[-1]
                
    subprocess.run('rm {0}'.format(fn), shell=True, check=True)
    
    return score
    
# Main code:
def main():
    """
    Main code.
    """
    # Define arguments:
    args = define_arguments().parse_args()
    fn = args.input
    
    # List GenBank IDs from Rosalind input file:
    with open(fn, 'r') as fo:
        ids = fo.read().split()
        
    # Retrieve records:
    #print('Retrieving records...')
    fasta_dict = {}
    for header, seq in retrieve_GenBank_record(ids):
        fasta_dict[header]= seq
        
    # Get two sequences from fasta_dict:
    keys = [key for key in fasta_dict.keys()]
    
    # Run alignment on sequences with Needle from EMBOSS:
    #print('Running alignment...')
    alignment_fn = run_needle(fasta_dict[keys[0]], fasta_dict[keys[1]])
    
    # Get alignment score from Needle output file:
    score = parse_needle(alignment_fn)
    print(score)
    
if __name__ == '__main__':
    main()
