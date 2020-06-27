#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving Rosalind exercise Majority Element
                    (maj).

Usage:              python3 maj.py <Rosalind input file>

Arguments:
    <Rosalidn input file>   .txt input file from Rosalind containing a
                            positive integer k and a positive integer n on the
                            first line, space separated. Consecutive lines
                            contain k arrays of n positive integers.
"""
from sys import argv

def parse_input(filename):
    """Parses two positive integers and several arrays from input file.

    Arguments:
    filename (str): name of Rosalind input file to parse from.

    Returns:
    k (int):        positive integer describing the number of arrays.
    n (int):        positive integer describing the length of all arrays.
    arrays (list):  list of lists of integers.

    Input file from Rosalind containing space separated positive integers n
    and k on the first line and each of the consecutive k lines an array
    of n integers.
    """
    arrays = []
    with open(filename, 'r') as file_open:
        k, n = list(map(int, file_open.readline().strip().split()))
        for line in file_open:
            arrays.append(list(map(int, line.strip().split())))
    return k, n, arrays

def majority_element(array):
    """Returns the majority element of an array, if it exists.

    Arguments:
    array (list):   list of integers to find majority element from.

    Returns:
    maj_elem (int): majority element from array, returns -1 if not exists.

    The input array should only contain positive numbers.
    """
    majority = int(len(array) / 2)
    elements = set(array)
    for element in elements:
        element_count = array.count(element)
        if element_count > majority:
            return element
    return -1

def main():
    # Parse input from Rosalind input file:
    k, n, arrays = parse_input(argv[1])
    # Find majority elements from input arrays:
    maj_elems = [majority_element(array) for array in arrays]
    # Print majority elements to stdout:
    print(" ".join(list(map(str, maj_elems))))

if __name__ == '__main__':
    main()
