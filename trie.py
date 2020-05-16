#!/usr/bin/env python3
"""
Author: David Meijer

Script for solving Rosalind assignment Introduction to Pattern Matching (trie).

Usage: "python3 trie.py <input_filename.txt>"
"""

# import statements
from sys import argv

# functions
def parse_input(filename):
    """Parses DNA sequences from input file and returns them in a list.
    
    filename: str, file to parse DNA sequences from.
    
    The input file contains a seperate DNA sequence on seperate line.
    """
    dnas = []
    
    with open(filename, 'r') as fo:
        for line in fo:
            dnas.append(line.strip())
        
    return dnas
    
class Node():
    """Class for creating a trie Node object.
    """
    def __init__(self, indx):
        """Initializes Node object.
        """
        self.indx = indx
        self.children = {}
        
    def __repr__(self):
        """Returns index string as Node representation.
        """
        return str(self.indx)
        
class Trie():
    """Class for creating a Trie object.
    """
    def __init__(self, origin = 1):
        """Initializes Trie object.
        """
        self.root = Node(origin)
        self.counter = origin
        
    def __repr__(self):
        """Returns string as Trie representation.
        """
        pass
        
    def insert(self, seq):
        """Inserts a new sequence into the Trie.
        """
        current_node = self.root
        
        for nucl in seq:
            if nucl not in current_node.children:
                self.counter += 1
                current_node.children[nucl] = Node(self.counter)
                
            current_node = current_node.children[nucl]
                    
    def edges(self, parent = False):
        """Recursive generator for edges from trie structure object.
        """
        if not parent:
            parent = self.root

        for edge, child in parent.children.items():
            yield parent, child, edge
            
            if child.children:
                for edge in self.edges(child):
                    yield edge
        
def make_trie(seqs):
    """Returns a Trie object created from input sequences.
    """
    trie = Trie()
    
    for seq in seqs:
        trie.insert(seq)
        
    return trie
    
# main
def main():
    
    # step 1: parse DNA sequences from input file 
    dnas = parse_input(argv[1])
    
    # step 2: make trie from DNA sequences
    trie = make_trie(dnas)
    
    # step 3: print edges from trie to stdout
    for parent, child, edge in trie.edges():
        print(parent, child, edge)
    
if __name__ == "__main__":
    main()
