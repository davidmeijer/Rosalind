#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving Rosalind exercise Reversal Distance
                    (rear).

Usage:              python3 rear.py <input.txt>

Arguments:
    input.txt       Rosalind assiggnment specific input file.
"""
from sys import argv

def parse_input(fi):
    """Generator for parsing permutation collections from input file."""
    with open(fi, "r") as fo:
        record = []
        for line in fo:
            if not line.strip():
                yield record
                record = []
                continue
            record.append(list(map(int, line.strip().split())))
        yield record

def get_lis(iter):
    """Returns longest increasing substring from iterable of integers."""
    n = len(iter)
    subs = [[] for _ in range(n)]
    subs[0].append(iter[0])
    for idx_a in range(1, n):
        for idx_b in range(idx_a):
            if iter[idx_a] > iter[idx_b] \
            and len(subs[idx_a]) < (len(subs[idx_b]) + 1):
                subs[idx_a] = subs[idx_b].copy()
        subs[idx_a].append(iter[idx_a])
    return max(subs, key = lambda sub: len(sub))

def reversal_dist(source, target):
    """Returns reversal distance between array source and array target."""
    source_map = {i:idx + 1 for (idx, i) in enumerate(source)}
    target = [source_map[j] for j in target]
    while source != target:


def main():
    """Driver code."""
    for pair in parse_input(argv[1]):
        count = reversal_dist(pair[0], pair[1])
        print(count)

if __name__ == "__main__":
    main()
