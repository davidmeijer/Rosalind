#!/usr/bin/env python3
"""
Author: David Meijer

Script to solve Rosalind exercise Fibonacci Numbers (fibo).

Run script in command line as: "python3 fibo.py <input_filename>"
"""

# import statements
from sys import argv

# functions
def parse_input(filename):
    """Returns parsed positive integer from input file.

    filename: str, name of input file.

    Input file contains one line with a single positive integer on the first
    line.
    """
    with open(filename, 'r') as fo:
        n = fo.readline().strip()

    return int(n)

def fibonacci(n):
    """Calculates up to n Fibonacci numbers.

    n: int, number of Fibonacci numbers to calculate.
    fibs: list of ints, list of n Fibonacci numbers that have been calculated.
    """
    if n == 0:
        return [0]

    fibs = [0, 1] + ([None] * (n - 1))

    for ind in range(2, len(fibs)):
        fibs[ind] = fibs[ind - 1] + fibs[ind - 2]

    return fibs

# main
def main():

    # step 1: parse input positive integer n from input file
    n = parse_input(argv[1])

    # step 2: return the nth Fibonacci number
    n_fibs = fibonacci(n)

    # step 3: print result to stdout
    print(n_fibs[n])

if __name__ == "__main__":
    main()
