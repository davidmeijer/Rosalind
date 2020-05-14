#!/usr/bin/env python3
"""
Author: David Meijer

Script for solving Rosalind exercise Counting Phylogentic Ancestors (inod).

Usage: "python3 inod.py <input_filename>"
"""

# import statements
from sys import argv

# functions
def parse_input(filename):
    """Parses a positive integer from input file.

    filename: str, name of input file containing a single positive integer.

    Input file contains a single positive integer on the first line.
    """
    with open(filename, 'r') as fo:
        n = fo.readline().strip()

    return int(n)

def internal_nodes(leaves):
    """Returns the number of internal leaves of unrooted tree with n leaves.
    
    leaves: int, number of leaves of unrooted tree.
    """
    return leaves - 2

# main
def main():

    # step 1: parse positive input from input file
    n = parse_input(argv[1])
    
    # step 2: count internal nodes
    nodes = internal_nodes(n)

    # step 3: print number of internal nodes to stdout
    print(nodes)

if __name__ == "__main__":
    main()
