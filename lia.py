#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse
from scipy.special import binom

def define_arguments():
    """Defines possible command line arguments.

    Returns:
        parser (object): contains command line input arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input file.')
    return parser

def parse_input(fn):
    """Parses Rosalind input file.

    Args:
        fn (str): Rosalind input file filename.

    Returns:
        input (list): list of two integers."""
    with open(fn) as fo:
        input = [int(x) for x in fo.readline().strip().split()]
    return input

def p(x, N):
    """Probability of at least N successes in x trials.

    Args:
        x (int): number of trials.
        N (int): number of successes.

    Returns:
        p_sum (float): chance of at least N successes."""
    p_sum = 1 - sum(Bernoulli(x, y) for y in range(N))
    return p_sum

def Bernoulli(x, y):
    """Calculates probability of y successes in x trials.

    Probability of AbAb from AbAb x AbAb is 4/16 or 0.25.
    Probability of non-AbAb is therefore 12/16 or 0.75.

    Args:
        x (int): number of trials.
        y (int): number of successes.

    Returns:
        p (float): chance of y successes in x trials."""
    p = binom(x, y) * (0.25 ** y) * (0.75 ** (x - y))
    return p

def main():
    args = define_arguments().parse_args()
    input = parse_input(args.input)
    result = round(p(2 ** input[0], input[1]), 3)
    print(result)

if __name__ == '__main__':
    main()
