#!/usr/bin/env python3
"""
Author: David Meijer

Script for solving Rosalind exercise Partial Permutations (pper).

More information on partial permutations: 
http://onlinestatbook.com/2/probability/permutations.html

Usage: "python3 pper.py <input_file.txt>"
"""

# import statements
from sys import argv

# functions
def parse_input(filename):
	"""Parses two positive integers from input file.
	
	filename: str, file name of input file containing two positive
	integers.
	
	Input file contains a single line with two positive integers
	seperated by a single space.
	"""
	with open(filename, 'r') as fo:
		n, k = fo.readline().strip().split()
		
	return int(n), int(k)
	
def factorial(n):
	"""Returns the factorial of n.
	
	n: int, positive integer.
	"""
	if n == 1:
		return 1
	
	else:
		return n * factorial(n - 1)
		
def nPr(n, r):
	"""Returns the number of combinations you can pick r from n.
	
	r: int, positive integer.
	n: int, positive integer.
	"""
	return int(factorial(n) / factorial(n - r))
	
# main
def main():
	
	# step 1: parse n and k from input file
	n, k = parse_input(argv[1])
	
	# step 2: calculate number of partial permutations, modulo 1,000,000
	pper = nPr(n, k) % 1000000

	# step 3: print result to stdout
	print(pper)
	
if __name__ == "__main__":
	main()
