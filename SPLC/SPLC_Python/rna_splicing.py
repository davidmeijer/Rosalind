#!/usr/bin/env python3
"""

Author: David Meijer

Rosalind exercise: RNA splicing

"""
import argparse
import string

codon_table = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    }


def define_arguments():
    """Defines possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parses exercise specific rosalind input.

    Args:
        path_to_input (str): absolute or relative path to file with
        specific input for this rosalind exercise.

    Returns:
        fasta_dict (dict): parsed sequences as {header : seq}.

    """
    fasta_dict = {}

    with open(path_to_input, 'r') as in_fo:
        for line in in_fo:
            line = line.strip()
            if len(line) > 0:
                if line.startswith('>'):
                    header = line[1:]
                    fasta_dict[header] = []
                    continue
                else:
                    fasta_dict[header].append(line)

    for header, seqs in fasta_dict.items():
        fasta_dict[header] = ''.join(seqs)

    return fasta_dict

def splicing(fasta_dict):
    """Splices longest seq in dict with smaller seqs.

    Args:
        fasta_dict (dict): parsed sequences as {header : seq}.

    Returns:
        spliced_dna (str): spliced DNA seq.

    """
    # Get longest seq from fasta_dict and accompanying header:
    longest_seq = ''
    longest_seq_header = ''
    for header, seq in fasta_dict.items():
        if len(seq) > len(longest_seq):
            longest_seq = seq
            longest_seq_header = header

    del fasta_dict[longest_seq_header]

    for header, seq in fasta_dict.items():
        partA, intron, partB = longest_seq.partition(seq)
        longest_seq = ''.join([partA, partB])

    return longest_seq

def transcribe(dna):
    """Transcribes DNA into RNA.

    Args:
        dna (str): DNA sequence.

    Returns:
        rna (str): transcribed DNA sequence into RNA.

    """
    trans_table = str.maketrans('GCTA', 'GCUA')

    rna = str.translate(dna, trans_table)

    return rna

def translate(rna):
    """Translates DNA into peptide sequence.

    Args:
        dna (str): spliced DNA seq.

    Returns:
        pep (str): translated spliced DNA seq.

    """
    pep = []

    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        pep.append(codon_table[codon])

    return ''.join(pep)


def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    fasta_dict = parse_rosalind_input(args.input)
    spliced_dna = splicing(fasta_dict)
    rna = transcribe(spliced_dna)
    pep = translate(rna)

    print(pep[:-1])

if __name__ == '__main__':
    main()