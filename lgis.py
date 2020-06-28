#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving Rosalind exercise Longest Increasing
                    Subsequence (lgis).

Usage:              python3 lgis.py <input.txt>

Arguments:
    <input.txt>     Rosalind specific input file.
"""
from sys import argv

def parse_input(fn):
    """Parses permutation and permutation length from input file."""
    with open(fn, "r") as fo:
        n = int(fo.readline().strip())
        iter = list(map(int, fo.readline().strip().split()))
    return n, iter

def get_lis(n, iter):
    """Returns longest increasing substring from iterable of integers."""
    subs = [[] for _ in range(n)]
    subs[0].append(iter[0])
    for idx_a in range(1, n):
        for idx_b in range(idx_a):
            if iter[idx_a] > iter[idx_b] \
            and len(subs[idx_a]) < (len(subs[idx_b]) + 1):
                subs[idx_a] = subs[idx_b].copy()
        subs[idx_a].append(iter[idx_a])
    return max(subs, key = lambda sub: len(sub))

def print_lis(lis):
    """Prints lis in nice format to stdout."""
    print(" ".join(list(map(str, lis))))

def main():
    """Driver code."""
    n, iter = parse_input(argv[1])
    print_lis(get_lis(n, iter))
    print_lis(get_lis(n, iter[::-1])[::-1])

if __name__ == "__main__":
	main()
