#!/usr/bin/env python3
"""Author: David Meijer

Description: Script takes two DNA strings s and t (each of length
at most 1 kbp) and returns all locations of t as a substring of s."""

#Space for imports:
from sys import argv
import re

#Space for definitions:
def parse_fo(fo):
	"""Parses input line 1 as DNA seq and line 2 as substring of DNA seq.

	fo -- file object, lines containing input fo.
	dna -- string, DNA sequence.
	sub --- string, DNA subsequence."""
	lines = []
	for line in fo:
		lines.append(line.strip())
	dna, sub = lines[0], lines[1]
	return dna, sub

def get_locs(dna, sub):
	"""Get start locations of sub in dna.

	dna -- string, DNA sequence.
	sub -- string, DNA subsequence.
	locs -- list of strings, start locations of sub in dna."""
	locs = []
	dna_cuts = []
	for i in range(len(dna)-1):
		dna_cuts.append(dna[i:i+len(sub)])
	results = dna_cuts.index(sub)
	for i, j in enumerate(dna_cuts):
		if j == sub:
			locs.append(str(i+1))
	return locs

#Space for main code:
if __name__ == "__main__":

	#Open fo and read line:
	with open(argv[1]) as fo:
		#Parse lines from fo:
		dna, sub = parse_fo(fo)

	#Get locations of sub as a substring of dna:
	locs = get_locs(dna, sub)

	#Open out file:
	with open(argv[1][:-4] + '_out.txt', 'w') as out_fo:
		#Print locations to out_fo:
		out_fo.write(' '.join(locs))
