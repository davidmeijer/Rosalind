#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Global Multiple Alignment
http://rosalind.info/problems/clus/
"""
# Imports:
import argparse
import subprocess
import copy

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

def parse_fasta(fn):
    """
    Parses FASTA input file into dict as {header:seq,...}.
    
    Args:
        fn (str): FASTA input file name.
    Returns:
        fasta_dict (dict): parsed FASTA input file as {header:seq,...}.
    """
    fasta_dict = {}
    
    with open(fn, 'r') as fo:
        for line in [line.strip() for line in fo]:
            if line.startswith('>'):
                header = line
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line
                
    return fasta_dict
    
def run_clustalw(fn):
    """
    Run multiple alignment on FASTA format sequences.
    
    Args:
        fn (str): file name of FASTA file with seqs to align.
    """
    cmd = ('clustalw {0} -PIM').format(fn, fn[:-4])
    subprocess.run(cmd, shell=True, check=True)
    
def hamm_clustalw_aln(aln_fn):
    """
    Calculates hamming distances between all aligned sequences.
    
    Args:
        aln_fn (str): file name of CLUSTALW output alignment file.
    Returns:
        header (str): header of sequence with highest combined hamming
                      score.
        i (int): combined hamming score.
    """
    fasta_dict = {}
    
    with open(aln_fn, 'r') as fo:
        for line in [line.strip() for line in fo][3:-1]:
            try:
                header, seq = line.split()
                if not header in fasta_dict:
                    fasta_dict[header] = seq
                else:
                    fasta_dict[header] += seq
            except:
                pass
    
    keys = [key for key in fasta_dict.keys()]
    for i in keys:
        hamm_score = 0
        for j in keys:
            if i != j:
                hamm_score += hamming(fasta_dict[i], fasta_dict[j])
        
        yield i, hamm_score

def hamming(seq1, seq2):
	"""
    Calculates the Hamming distance between two sequences.
	
    Args:
        seq1 (str): DNA sequence 1.
        seq2 (str): DNA sequence 2.
        
    Returns:
        hamm (int): hamming distance between seq1 and seq2.
	"""
	hamm = 0
	for i in range(len(seq1)-1):
		if seq1[i] != seq2[i]:
			hamm += 1

	return hamm
    
# Main code:
def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    fasta_dict = parse_fasta(args.input)
    
    # Parse input file so that sequence is one one line for CLUSTALW:
    fi = '../input/parse_fasta.txt'
    with open(fi, 'w') as fi_fo:
        for header, seq in fasta_dict.items():
            fi_fo.write('{0}\n{1}\n'.format(header, seq))
            
    run_clustalw(fi)
    
    parse_fasta(fi)
    outfn = fi[:-4] + '.aln'
    
    for header, score in hamm_clustalw_aln(outfn):
        print(header, score)
    
    
if __name__ == '__main__':
    main()
