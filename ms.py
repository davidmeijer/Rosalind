#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving the Rosalind assignment Merge Sort (ms).

Usage:              python3 ms.py <Rosalind input file>

Arguments:
    <Rosalind input file>   .txt input file from Rosalind containing a
                            positive integer on the first line and an
                            unsorted array of integers on the second line.
"""
from sys import argv

def parse_input(filename):
    """Parses a positive integer and an unsorted list from input file.

    Arguments:
    filename (str): Rosalind assignment input file to parse input from.

    Returns:
    n (int):        length of the unsorted array of integers array_n.
    array_n (list): unsorted array of integers with length n.

    Rosalind input file contains a positive integer on the first line and an
    unsorted array of integers on the second line.
    """
    with open(filename, 'r') as file_open:
        n = int(file_open.readline().strip())
        array_n = list(map(int, file_open.readline().strip().split()))

    return n, array_n

def merge_sort(array_n):
    """Returns a merge sorted array.

    Arguments:
    array_n (list): unsorted array of integers with length n.

    Returns:
    sorted (list): returns the sorted array of integers with length n.
    """
    if len(array_n) <= 1:
        return array_n

    left, right = [], []
    for index, number in enumerate(array_n):
        if index < (len(array_n) / 2):
            left.append(number)
        else:
            right.append(number)

    return merge(merge_sort(left), merge_sort(right))

def merge(left, right):
    """Returns a merged list.

    Arguments:
    left (list):    sorted list of integers.
    right (list):   sorted list of integers.

    Returns:
    result (list):  sorted list of integers, merged left and right.
    """
    result = []

    while len(left) > 0 and len(right) > 0:

        if left[0] <= right[0]:
            j = left.pop(0)
            result.append(j)
        else:
            j = right.pop(0)
            result.append(j)

    if len(left) > 0:
        result += left

    if len(right) > 0:
        result += right

    return result

def main():
    # Parse input from Rosalind input file:
    n, array_n = parse_input(argv[1])

    # Merge sort the unsorted array_n with length n:
    sorted_array = merge_sort(array_n)

    # Print sorted array to stdout:
    print(" ".join(list(map(str, sorted_array))))

if __name__ == '__main__':
    main()
