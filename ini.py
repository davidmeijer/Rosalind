#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
from Bio.Seq import Seq

def define_arguments():
	"""Defines possible command line arguments.

	Returns:
		parser (obj): contains user input command line arguments."""
	parser = argparse.ArgumentParser()
	parser.add_argument('-input', type=str, required=True,
						help='Rosalind input file.')
	
	return parser

def main():
	"""Main code."""
	args = define_arguments().parse_args()
	
	with open(args.input, 'r') as fo:
		seq = fo.read().strip()

	print(seq.count('A'), 
          seq.count('C'), 
          seq.count('G'), 
          seq.count('T'))

if __name__ == '__main__':
	main()
