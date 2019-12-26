#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: LCSQ
http://rosalind.info/problems/lcsq/
"""
import argparse
import itertools
import random
import pandas
import sys

sys.setrecursionlimit(1500)

def define_arguments():
    """
    Defines possible command line arguments.

    Returns:
        parser (object): contains user arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
    return parser

def parse_rosalind_input(fn):
    """
    Parses Rosalind input.

    Args:
        fn (str): file name of Rosalind input file.
    """
    fasta_dict = {}
    with open(fn) as fo:
        for line in [line.strip() for line in fo]:
            if line.startswith('>'):
                header = line
                fasta_dict[header] = []
            else:
                fasta_dict[header].append(line)
    fasta_dict.update(
        (header, ''.join(seq)) for header, seq in fasta_dict.items())
    return fasta_dict

def get_subseqs(dna):
    """
    Creates different ascending subseqs from DNA seq.

    Args:
        dna (str): DNA sequence.

    Returns:
        subseqs (list): list of sub DNA seqs from input seq.
    """
    subseqs = [dna] # For if the two seqs are identical they match in full length.
    n = [i+1 for i in range(len(dna))]
    for subset in range(len(dna)):
        for lst in [list(j) for j in
                    [x for x in itertools.permutations(n, subset+1)]]:
            if order(lst):
                # Now you know which bases to omit, construct subseqs:
                subseq = []
                for i,base in enumerate(dna):
                    if i+1 not in lst:
                        subseq.append(base)
                subseqs.append(''.join(subseq))
    return list(set(subseqs))

def randomString(stringLength):
    """Generate a random string of fixed length """
    bases = 'ACTG'
    return ''.join(random.choice(bases) for i in range(stringLength))

def traceback(df, loc, route=[]):
    """"""

    """
    x = df.at[loc[0]-1,loc[1]-1]
    y = df.at[loc[0]-1,loc[1]]
    z = df.at[loc[0],loc[1]-1]
    c = df.at[loc[0], loc[1]]

    #print(loc, x, y, z, route)

    #if loc == (1,1):
    if x == 0 and y == 0 and z == 0:
        route.append(loc[1])
        return route

    if x == y and x == z and x == c:
        return traceback(df, (loc[0]-1,loc[1]-1), route)

    if x == y and x == z and c > x:
        route.append(loc[1])
        return traceback(df, (loc[0]-1,loc[1]-1), route)

    if y > z:
        return traceback(df, (loc[0]-1,loc[1]), route)

    if z > y:
        return traceback(df, (loc[0],loc[1]-1), route)

    # At intersections you can go in two ways, I just picked one which
    # was not the longest...
    #if z == y and x < z:
    #    return traceback(df, (loc[0],loc[1]-1), route)

    if z == y and x < z:
        return traceback(df, (loc[0]-1,loc[1]), route)
    """

def main():
    """
    Main code.
    """
    args = define_arguments().parse_args()
    fasta_input = parse_rosalind_input(args.input)

    seqs = [fasta_input.get([*fasta_input][i]) for i in range(len(fasta_input))]

    #s = seqs[0]
    #t = seqs[1]

    #s = 'AGCAT'
    #t = 'GAC'

    s = 'TTTTTAAAAA'
    t = 'AAAAATTTTT'

    #s = randomString(1000)
    #t = randomString(1000)

    #print(len(s), len(t))

    df = pandas.DataFrame(columns=[0] + [x+1 for x in range(len(s))],
                          index=[0] + [x+1 for x in range(len(t))])
    df.iloc[0] = [0]; df[0]= 0

    for i,row in enumerate(df.index):
        if row == 0:
            pass
        else:
            for j,column in enumerate(df):
                if column == 0:
                    pass
                else:
                    x = max([df.at[i,j-1],df.at[i-1,j]])
                    if df.at[i,j-1] == df.at[i-1,j]:
                        if s[j-1] == t[i-1]:
                            x += 1
                    df.at[i,j] = x

    print(df)

    start = (df.shape[0]-1, df.shape[1]-1)
    route = traceback(df, start)

    print(route)
    #print(df)
    lcs = ''.join([s[loc-1] for loc in route][::-1])
    print(lcs)
    print(len(lcs))

if __name__ == '__main__':
    main()