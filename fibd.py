#!/usr/bin/env python3
"""
Author: David Meijer

Script for solving Rosalind exercise Mortal Fibonacci Rabbits (fibd).

Usage: "python3 fibd.py <input_filename.txt>"
"""

# import statements
from sys import argv

# functions
def parse_input(filename):
	"""Parses positive integers n and m from input file.
	
	filename: str, input file containing positive integers n and m.
	n: int, number of months rabbits are counted.
	m: int, number of months all rabbits live.
	
	Input file contains a single line with two positive integers as n m. 
	"""
	with open(filename, 'r') as fo:
		n, m = fo.readline().strip().split()
		
	return int(n), int(m)
	
def mortal_fibonacci(months, life_exp):
	"""Returns the number rabbit pairs alive after the n-th month.
	
	months: int, number of month cycles rabbits are counted.
	live_exp: int, number of months all rabbits live.
	"""
	mature = [0] * life_exp
	infant = [0] * life_exp
	
	infant[-1] = 1
	
	for cycle in range(0, months -1):
		mature, infant = add_generation(mature, infant)
	
	return mature [-1] + infant[-1]
		
def add_generation(mature, infant):
	"""Adds generation to mature and infant lines.
	"""
	mature.append(mature[-1] + infant[-1] - infant[0])
	infant.append(mature[-2])
	
	return mature[1:], infant[1:]

# main
def main():
	
	# step 1: parse positive integers n and m from input file
	n, m = parse_input(argv[1])
	
	# step 2: calculate the number of rabbits that will remain after the
	#		  n-th month if all rabbits live for m months
	number_of_rabbits = mortal_fibonacci(n, m)
	
	# step 3: print the resulting number of rabbits to stdout
	print(number_of_rabbits)
	
if __name__ == "__main__":
	main()
