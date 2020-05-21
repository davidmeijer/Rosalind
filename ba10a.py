#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Script for solving Rosalind exercise Compute the
                    Probability of a Hidden Path (ba10a).

Usage:              python3 ba10a.py <Rosalind input file>

Arguments:
<Rosalind input file>   .txt input file containing a hidden path on the first
                        line, the states on the third line and from line
                        three on the transition matrix.
"""
from sys import argv

def parse_input(filename):
    """Script for parsing the Rosalind assignment input data.

    Arguments:
    filename (str):     name of Rosalind assignment input file.

    Returns:
    hidden_path (str):  the hidden path parsed from input file.
    states (list):      list of strings describing the states in the path.
    matrix (list):      list of list describing the transition HMM matrix.
    """
    matrix = []

    with open(filename, 'r') as fo:
        sep_count = 0

        for line in fo:
            if line.startswith('---'):
                sep_count += 1
                continue

            if sep_count == 0:
                hidden_path = line.strip()

            if sep_count == 1:
                states = line.strip().split()

            if sep_count == 2:
                matrix_line = line.strip().split()
                for idx, item in enumerate(matrix_line):
                    try:
                        item = float(item)
                        matrix_line[idx] = item
                    except:
                        pass
                matrix.append(matrix_line)

    return hidden_path, states, matrix

def probability_path(path, matrix):
    """Returns all probabilities to go from node to node in hidden path.

    Arguments:
    path (str):         hidden path.
    matrix (list):      list of lists describing the transition probs.

    Returns:
    probs (list):       returns transition probs for hidden path.
    """
    probs = [1 / len(matrix[0])]

    for i in range(len(path)):
        if i != (len(path) - 1):
            start = matrix[0].index(path[i])
            final = matrix[0].index(path[i + 1])
            prob = matrix[start + 1][final + 1]
            probs.append(prob)

    return probs

def multiply_list(array):
    """Multiplies all numbers in list.

    Arguments:
    array (list):       list of integers to multiply.

    Returns:
    result (int):       multiplication of all integers in array.
    """
    while len(array) > 1:
        a, b = array.pop(0), array.pop(0)
        array.append(a * b)

    result = array[0]

    return result

def main():
    # Parse input from Rosalind file:
    path, states, matrix = parse_input(argv[1])

    # Calculate the probability of the path:
    probs = probability_path(path, matrix)
    prob = multiply_list(probs)

    # Print calculated probability to stdout:
    print(prob)

if __name__ == '__main__':
    main()
