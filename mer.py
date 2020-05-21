#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving the Rosalind assignment Merge Two
                    Sorted Arrays (mer).

Usage:              python3 mer.py <Rosalind input file>

Arguments:
    <Rosalind input file>   .txt input file from the Rosalind assignment
                            containing a positive integer n on the first
                            line, an array of length n on the second line,
                            a positive integer m on the third line, and an
                            array of length m on the fourth line.
"""
from sys import argv
import sys

def parse_input(filename):
    """Input file containing two sorted arrays and their lengths.

    Arguments:
    filename (str): Rosalind assignment input file to parse input from.

    Returns:
    n (int):        positive integer describing length of sorted array_n.
    m (int):        positive integer describing length of sorted array_m.
    array_n (list): list with length n containing sorted integers.
    array_m (list): list with length m containing sorted integers.

    Input is from Rosalind containing a positive integer n on the first line,
    an array of length n on the second line, a positive integer m on the third
    line, and an array of length m on the fourth line.
    """
    with open(filename, 'r') as file_open:
        n = int(file_open.readline().strip())
        array_n = list(map(int, file_open.readline().strip().split()))
        m = int(file_open.readline().strip())
        array_m = list(map(int, file_open.readline().strip().split()))

    return n, m, array_n, array_m

def merge_sorted_arrays(array_n, array_m, merged=[]):
    """Returns a sorted array constructed from two sorted arrays.

    Arguments:
    array_n (list): list with length n containing sorted integers.
    array_m (list): list with length m containing sorted integers.

    Returns:
    merged (list):  merged list with length n + m.
    """
    if len(merged) == 0:
        array_n.append('*')
        array_m.append('*')

    if array_n[0] != '*' or array_m[0] != '*':
        if array_n[0] <= array_m[0]:
            merged.append(array_n.pop(0))
        else:
            merged.append(array_m.pop(0))

    if array_n[0] == '*':
        merged += array_m[:-1]
        return merged

    if array_m[0] == '*':
        merged += array_n[:-1]
        return merged

    return merge_sorted_arrays(array_n, array_m, merged)

def main():
    sys.setrecursionlimit(20000)

    # Parse input from Rosalind exercise specific input file:
    n, m, array_n, array_m = parse_input(argv[1])

    # Merge the two parsed sorted arrays:
    merged_arrays = merge_sorted_arrays(array_n, array_m)

    # Print merged sorted array to stdout:
    print(" ".join(list(map(str, merged_arrays))))

if __name__ == '__main__':
    main()
