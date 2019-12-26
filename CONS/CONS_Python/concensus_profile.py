#!/usr/bin/env python3
"""

Author: David Meijer

Rosalind exercise: Concensus and Profile

"""
import argparse
import pandas as pd

def define_arguments():
    """Define possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parses specific rosalind input file.

    Args:
        path_to_input (str): relative or absolute path to input file
        containing rosalind input data.

    Returns:
        fasta_dict (dict): dictionary containing input DNA fasta
        sequences.

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

def create_df(fasta_dict):
    """Creates dataframe for storing ACTG count.

    Args:
        fasta_dict (dict): dictionary containing input DNA fasta
        sequences.

    Returns:
        df (pandas dataframe): contains ACTG count per position.

    """
    # Create four rows (ACTG) and approriate number of columns per
    # sequence string position:
    index = ['A', 'C', 'T', 'G']
    for header, seq in fasta_dict.items():
        columns = []
        for i in range(len(seq)):
            columns.append(i)
        break

    df = pd.DataFrame(0, index, columns)

    for pos, column in df.iteritems():
        for header, seq in fasta_dict.items():
            actg_id = seq[pos]
            df.at[actg_id, pos] += 1

    return df

def create_concensus_seq(df):
    """Creates concensus seq based on ACTG counts.

    Args:
        df (pandas dataframe): contains ACTG count per position.

    Returns:
        concensus_seq (str): concensus sequence based on counts.

    """
    concensus_seq = []

    for max_count  in df.idxmax(axis=0):
        concensus_seq.append(max_count)

    return ''.join(concensus_seq)


def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    fasta_dict = parse_rosalind_input(args.input)
    df = create_df(fasta_dict)
    concensus_string = create_concensus_seq(df)

    print(concensus_string)
    print('A: ' + ' '.join([str(x) for x in df.loc['A']]))
    print('C: ' + ' '.join([str(x) for x in df.loc['C']]))
    print('G: ' + ' '.join([str(x) for x in df.loc['G']]))
    print('T: ' + ' '.join([str(x) for x in df.loc['T']]))

if __name__ == '__main__':
    main()