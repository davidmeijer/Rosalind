#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
import subprocess

def define_arguments():
    """
    Defines possible command line arguments.
    
    Returns:
        parser (obj): contains command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
                        
    return parser
    
def fastq_to_fasta(fn_in, 
                   fn_out=None, 
                   del_old=False):
    """Converts FASTQ file to FASTA format.
    
    Args:
        fn_in (str): input file name for file in FASTQ format.
        fn_out (str): output file name for file in FASTA format.
        del_old (bool): delete input file."""
    if fn_out == None:
        # Only split on last '.' delimiter to exchange format:
        fn_out = ''.join(fn_in.rsplit('.', 1)[:-1]) + '.fq'
        
    with open(fn_in, 'r') as in_o, open(fn_out, 'w') as out_o:
        entry_count = 0 # Every entry contains four lines in FASTQ.
        for line in [line.strip() for line in in_o]:
            entry_count += 1
            if entry_count == 1:
                header = '>' + line[1:]
                out_o.write(header + '\n')
            if entry_count == 2:
                sequence = line
                out_o.write(sequence + '\n')
            if entry_count == 3:
                pass
            if entry_count == 4:
                entry_count = 0
                
    if del_old:
        cmd = 'rm {0}'.format(fn_in)
        subprocess.run(cmd, shell=True, check=True)

def main():
    args = define_arguments().parse_args() 
    fastq_to_fasta(args.input)

if __name__ == '__main__':
    main()
