#!/usr/bin/env python3
"""
Author: David Meijer
Assignment: INtroduction to Protein Databases
http://rosalind.info/problems/dbpr/
"""
# Imports:
import argparse
import urllib.request
import time
import re

#from Bio import ExPASy
#from Bio import SwissProt

# Classes and functions:
def define_arguments():
    """
    Defines possible command line arguments.
    
    Returns:
        parser (obj): contains user input command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
                        
    return parser
    
class UniProt:
    """
    For requesting UniProt entries.
    """
    def __init__(self, uniprot_id=None):
        self.up_id = uniprot_id
        
    def request_entry(self,
                      up_url="http://www.uniprot.org/uniprot/",
                      up_format=".txt",
                      up_random="?query=&random=yes"):
        """
        Request entry from UniProt with UniProt identifier.
        
        Dependencies: urllib.request, time, re
        
        Args:
            up_id (str): UniProt identifier.
            up_url (str): UniProt web address.
            up_format (str): retrieved format.
            up_random (str): if up_id=None, retrieve random entry.
            
            Available formats: txt,xml,rdf,fasta,gff
        
        Returns:
            raw_data (str): retrieved raw UniProt entry.
            up_id (str): UniProt identifier.
        """
        self.format = up_format
        
        try:
            if self.up_id != None:
                url = up_url + self.up_id + up_format
                raw_data = urllib.request.urlopen(url).read().decode("utf-8")
                self.rd = raw_data
            else:
                url = up_url + up_random
                html = urllib.request.urlopen(url).read().decode("utf-8")
                regex = re.compile("var entryId = '(.+)'")
                self.up_id = regex.search(html).group(1)
                time.sleep(.1)
                self.request_entry()
        except:
            print("Could not retrieve UniProt entry!")
            exit()
            
    def GO(self):
        """
        Retrieve Gene Ontology from raw UniProt entry.
        
        Args:
            raw_data (str): retrieved raw UniProt entry.
            up_format (str): retrieved format.
            
        Returns:
            GO (list): Gene Ontology terms as strings in list.
        """
        GO = []
        
        if self.format == '.txt':
            data = self.rd.split('\n')
            del self.rd, self.format
            for line in [line.split() for line in data]:
                if 'DR' in line and 'GO;' in line:
                    GO.append(line[4])
                    
        self.GO = GO

def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    
    with open(args.input, 'r') as fo:   
        UniProt_id = fo.read().strip()
        
    #handle = ExPASy.get_sprot_raw(UniProt_id)
    #record = SwissProt.read(handle)
    #dir(record)
    
    protein = UniProt(UniProt_id)
    protein.request_entry()
    protein.GO()
    
    bases = []
    branches = []
    for term in protein.GO:
        if not ';' in term:
            bases.append(term)
        else:
            branches.append(term[:-1])
    
    if len(bases) == 0:
        for j in range(len(branches)):
            print(branches[j])
    else:
        for i in range(len(bases)):
            for j in range(len(branches)):
                print(bases[i], branches[j])

if __name__ == '__main__':
	main()
