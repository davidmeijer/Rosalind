#!/usr/bin/env python3

"""Name: David Meijer

Usage: python3 rosalind_grph.py rosalind_grph.txt

Description: takes a collection of DNA strings in FASTA format having
a total length at most 10 kbp and returns the adjaceny list
corresponding to O(3). The edges are returned in random order."""

#Space for modules:
from sys import argv

#Space for definitions:
def parsing_fasta(infile):
    """Parses DNA strings in FASTA format into dictionary.

    infile -- file object, lines describing multiplte DNA strings in
    FASTA format."""
    seq_dict = {}
    for line in infile:
        if line.startswith(">"):
            fasta_id = line[1:].strip()
            seq_dict[fasta_id] = ""
        else:
            seq_dict[fasta_id] += line.strip()
    return seq_dict

def trim_seqs(seq_dict, adjacency_length=3):
    """Makes dictionary of FASTA seqs as {seq_id, trim_start, trim_end}.

    seq_dict -- dictionary, as {seq_id, seq}.
    seq_dict_out -- dictionary, as {seq_id, trim_start, trim_end}."""
    for fasta_id, seq in seq_dict.items():
        trim_seq_beg = seq[0:adjacency_length]
        trim_seq_end = seq[-adjacency_length:]
        seq_dict[fasta_id] = [trim_seq_beg,trim_seq_end]
    seq_dict_out = seq_dict
    return seq_dict_out

def finding_adjacency(seq_dict):
    """Finds adjacency between two sequences.

    seq_dict -- dictonary, dictionary, as {seq_id, trim_start, trim_end}.
    result_list -- list, with adjacencies as [seq_id1, seq_id2]."""
    result_list = []
    compare_dict = seq_dict.copy()
    for fasta_id, trim_list in seq_dict.items():
        check = [fasta_id, trim_list[0], trim_list[1]]
        for compare_id, compare_trim_list in compare_dict.items():
            if check[1] == compare_trim_list[1] and compare_id != check[0]:
                if not [compare_id, check[0]] in result_list:
                    result_list.append([compare_id, check[0]])
            if check[2] == compare_trim_list[0] and compare_id != check[0]:
                if not [check[0], compare_id] in result_list:
                    result_list.append([check[0], compare_id])
    return result_list

#Space for main code:
if __name__ == "__main__":

    #Step 1: open infile:
    with open(argv[1]) as fo:

        #Step 2: parse infile:
        fasta_dict = parsing_fasta(fo)

        #Step 3: trim fasta_dict:
        trim_dict = trim_seqs(fasta_dict)

        #Step 4: find adjacency:
        adjacencies = finding_adjacency(trim_dict)

    #Step 5: print results:
    with open(argv[1][:-4] + '_out.txt', 'w') as out_fo:
    	for result in adjacencies:
    		out_fo.write("{0} {1}\n".format(result[0], result[1]))
