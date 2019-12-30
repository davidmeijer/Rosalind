#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: Suboptimal Local Alignment
http://rosalind.info/problems/subo/
"""
# Imports:
import argparse
import subprocess


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
    Parses input fasta file in dictionary with {header:seq,...}.
    
    Args:
        fn (str): name of input file in FASTA format.
    
    Returns:
        fasta_dict (dict): parsed FASTA file as {header:seq,...}.
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
    
def run_lalign(fn1, fn2, outfn='out.txt'):
    """
    Run LALIGN on two sequences to find inexact repeat r (100% match).
    
    Args:
        fn1 (str): FASTA file with one sequence.
        fn2 (str): FASTA file with one sequence.
        
    Returns:
        outfn (str): file name of output file of LALIGN.
    """
    cmd = ('lalign36 -d 1 -f 8 -g 8 -m 3 ' +
           '-O {0} {1} {2}').format(outfn, fn1, fn2)
    subprocess.run(cmd, shell=True, check=True)
    
    return outfn
    
def parse_lalign(fn):
    """
    Retrieves 100.0% match from LALIGN alignment output file.
    
    Args:
        fn (str): file name LALIGN alignment output file.
        
    Returns:
        r (str): inexact repeat from input sequences.
    """
    r = ''
    
    line_count = 0
    with open(fn, 'r') as fo:
        for line in [line.strip() for line in fo]:
            line_count += 1
            if line.startswith('100.0% identity'):
                reset = line_count
            try:
                if line_count == reset + 2:
                    return line
            except:
                pass
                
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
    
def count_hamming(r, seq, max_hamm=3):
    """
    Count number of matches within dist mismatches.
    
    Args:
        r (str): inexact repeat from input sequences.
        seq (str): sequence to match r against.
        diost (int):
    Returns:
        count (int):
    """
    count = 0
    for i in range(len(seq) - len(r) + 1):
        if hamming(seq[i:i + len(r)], r) <= max_hamm:
            count += 1
            
    return count
    
# Main code:
def main():
    """
    Main code.
    """
    # Define input:
    args = define_arguments().parse_args()
    fn = args.input
    
    # Parse input file:
    fasta_dict = parse_fasta(fn)
    
    fns = []
    # Create seperate input files for input sequences:
    for header, seq in fasta_dict.items():
        fn = header[1:] + '.txt'
        fns.append(fn)
        with open(fn, 'w') as fo:
            fo.write('{0}\n{1}\n'.format(header, seq))
    
    # Find r:
    outfn = run_lalign(fns[0], fns[1])
    r = parse_lalign(outfn)
    
    # Count matches (max hamming distance is 3):
    for header, seq in fasta_dict.items():
        count = count_hamming(r, seq, 3)
        print(header, count)
    
if __name__ == '__main__':
    main()
