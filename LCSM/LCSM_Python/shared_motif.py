#!/usr/bin/env python3
"""

Author: David Meijer

Rosalind exercise: Find a Shared Motif.

"""
import argparse
import random
import copy

def define_arguments():
    """Defines possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parse specific rosalind input for this exercise.

    Args:
        path_to_input (str): absolute or relative path to input file
        with Rosalind exercise specific input lines in fasta format.

    Returns:
        fasta_dict (dict): containing parsed input fasta lines as
        {header : seq, ...}.

    """
    fasta_dict = {}

    with open(path_to_input, 'r') as in_fo:
        for line in in_fo:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                fasta_dict[header] = []
                continue
            else:
                fasta_dict[header].append(line)

    for header, seqs in fasta_dict.items():
        fasta_dict[header] = ''.join(seqs)

    return fasta_dict


class find_shared_motif():
    """Functions for finding largest shared motif among all fasta seqs.

    """
    def __init__(self, fasta_dict):
        """Class input arguments.

        Args:
            fasta_dict (dict): containing parsed input fasta lines as
            {header : seq, ...}.

        Returns:
            shared_motifs (list): all shared motifs found.

        """
        self.fasta_dict = fasta_dict
        self.shared_motifs = []


    def pipeline(self):
        """Pipeline for finding largest shared motif.

        """
        # Find largest shared motif between two random seqs and trim of
        # other part of seqs:
        header1, seq1 = self.choose_random_entries()
        header2, seq2 = self.choose_random_entries()
        motifs = self.find_shared_motif(seq1, seq2)

        omnipresent_motifs = []

        for motif in motifs:
            motif_count = 0
            for header, seq in self.fasta_dict.items():
                if motif in seq:
                    motif_count += 1
            if motif_count == len(self.fasta_dict):
                omnipresent_motifs.append(motif)

        print(max(omnipresent_motifs, key=len))

        # Determine largest shared motif:


    def choose_random_entries(self):
        """Chooses entry from dictionary at random and deletes it.

        Returns:
            header (str): header of random DNA sequence.
            seq (str): random DNA sequence.

        """
        header, seq = random.choice(list(self.fasta_dict.items()))
        del self.fasta_dict[header]

        return header, seq


    def find_shared_motif(self, seq1, seq2):
        """Find largest shared motif between two seqs.

        Args:
            seq1 (str): DNA sequence 1.
            seq2 (str): DNA sequence 2.

        Returns:
            motif (str): largest shared motif between seq1 and seq2.

        """
        slices = self.partition(seq1)

        slice_matches = []

        for slice in slices:
            if slice in seq2:
                slice_matches.append(slice)

        slice_matches = list(set(slice_matches))

        return slice_matches

    def partition(self, s):
        """Partitions string into all possible substrings.

        Args:
            s (str): string of any length.

        Returns:
            substrings (list): all unique possible substrings of s.

        """
        substrings = []

        for j in range(len(s)):
            ss = s[j:]
            for i in range(len(ss)):
                front, back = ss[:i], ss[i:]
                substrings.append(front)
                substrings.append(back)

        return list(set(substrings))

    def check_entries_on_motif(self, motif):
        """Delete entries with motif match.

        Args:
            motif (str): largest shared motif until now.

        Returns:
            motif (str): shorter largest shared motif..

        """
        pop_list = []

        for header, seq in self.fasta_dict.items():
            if seq.find(motif) == 0:
                pop_list.append(header)

        for header in pop_list:
            del self.fasta_dict[header]

def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    fasta_dict = parse_rosalind_input(args.input)
    find_shared_motif(fasta_dict).pipeline()


if __name__ == '__main__':
    main()