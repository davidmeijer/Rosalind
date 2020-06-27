#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving Rosalind exercise Binary Search (bins).

Usage:              python3 bins.py <Rosalind input file>

Arguments:
    <Rosalind input file>       .txt input file containing an integer n on
                                the first line, an integer k on the second
                                line, an array of integers of length n on the
                                third line, and an array of integers of length
                                k on the fourth line.
"""
from sys import argv

def parse_input(filename):
    """Parses two input arrays and their lengths from the input file.

    Arguments:
    filename (str): name of Rosalind input file.

    Returns:
    n (int):        positive integer describing length of array n.
    k (int):        positive integer describing length of array k.
    array_n (list): list of integers with length n.
    array_k (list): list of integers with length k.

    Rosalind exercise specific input file has n on the first line, k on the
    second line, array_n on the third line, and array_k on the fourth line.
    """
    with open(filename, 'r') as file_open:
        n = int(file_open.readline().strip())
        k = int(file_open.readline().strip())
        array_n = list(map(int, file_open.readline().strip().split()))
        array_k = list(map(int, file_open.readline().strip().split()))
    return n, k, array_n, array_k

def binary_search(target, array, index=1):
    """Returns one-based index of target in array.

    Arguments:
    target (int):   positive integer describing the target to find index for.
    array (list):   list of integers to find the target in.
    index (int):    found index of target in array (default: 1).

    Returns:
    index (int):    denoting the one-based index of target in array.

    The function will return the integer -1 if target is not found in array.
    """
    if len(array) > 1:
        split = int(len(array) / 2)
        left, right = array[:split], array[split:]
        if target in left:
            return binary_search(target, left, index)
        else:
            return binary_search(target, right, index + len(left))
    if len(array) == 1:
        if target == array[0]:
            return index
        else:
            return -1

def main():
    # Parse n, k, array_n, and array_k from input file:
    n, k, array_n, array_k = parse_input(argv[1])
    # Perform a binary search for all numbers in array_k:
    inds = [binary_search(i, array_n) for i in array_k]
    # Print indices to stdout:
    print(" ".join(list(map(str, inds))))

if __name__ == '__main__':
    main()
