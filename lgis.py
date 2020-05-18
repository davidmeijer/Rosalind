#!/usr/bin/env python3
"""
Author:         David Meijer
Description:    Script for solving Rosalind exercise Longest Increasing 
                Subsequence (lgis).

Usage:          python3 lgis.py <Rosalind input file>

Arguments:
    <Rosalind input file>   .txt input file containing a single positive 
                            integer on the first line and a collection of
                            space separated positive integers on the second 
                            line.
"""
# Import statements
from sys import argv

# Functions
def parse_input(filename):
	"""Parses positive integer and permutation from input file.
	
    Arguments;
        filename: str, file containing input integer and permutation.
        
    Returns:
        permutation_length: int, length of the permutation.
        permutation: list of ints, describing a permutation of length 
                     permutation_length.
	
	Input file consists of two lines: on the first line there is a 
	single integer and on the second line there is a permutation, 
	several numbers seperated each by a space.
	"""
	with open(filename, "r") as file_open:
		permutation_length = int(file_open.readline().strip())
		permutation = list(map(int, file_open.readline().strip().split()))
		
	return permutation_length, permutation
    
def create_graph(permutation):
    """Returns directed graph from permutation.
    
    Arguments:
        permutation_length: int, length of the permutation.
        permutation: list of ints, describing a permutation of length 
                     permutation_length.  
                     
    Returns:
        graph: dict, describing directions from number to number (nodes)
               that can be taken to get an ascending path through the
               permutation as {parent node (int) : [child nodes (int], ...}.
    """            
    graph = {}
    for value in permutation:
        graph[value] = []
        
    for parent in graph.keys():
        for value in permutation:
            if permutation.index(value) > permutation.index(parent) \
            and value > parent:
                graph[parent].append(value)
                
    return graph

def longest_path(graph):
    """Finds longest ascending path in graph.
    
    Arguments:
        graph: dict, describing directions from number to number (nodes)
               that can be taken to get an ascending path through the
               permutation as {parent node (int) : [child nodes (int], ...}.
               
    Returns:
    
    """
    # keep some kind of queue when walking through graph, when node has
    # no more neighbours and unvisited nodes in queue is also empty then
    # return path ...
    paths = []
    
    

# Main
def main():
    
    # Define name Rosalind input file
    input_filename = argv[1]
	
	# Step 1: parse permutation length and permutation from input file
    permutation_length, permutation = parse_input(input_filename)

	# Step 2: create graph from the permutaion for ascending routes
    graph = create_graph(permutation)
    
    # Step 3: find longest ascending path in graph
    path = longest_path(graph)
    
    # Step 4: print longest ascending path in permutation to stdout
    print(path)
    
if __name__ == "__main__":
	main()
