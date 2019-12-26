#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: Enumerating Oriented Gene Orderings
"""
import argparse
import itertools

def define_arguments():
    """Define possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Path to Rosalind exercise specific ' +\
                        'input file.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parses Rosalind exercise specific input file lines.

    Args:
        path_to_input (str): path to Rosalind exercise specific input
        file.

    Returns:
        pos_integer (int): positive integer n <= 6.

    """
    with open(path_to_input, 'r') as in_fo:
        pos_integer = int(in_fo.readline().strip())

    return pos_integer

def get_signed_permutations(pos_integer):
    """Get all signed parmutations from positive integer.

    Args:
        pos_integer (int): positive integer n <= 6.

    Returns:
        signed_perms (list): list of list as [ [signed_perm], ...].

    """
    signed_perms = []

    # Get all present numbers part of the signed permutations:
    elements = []
    for i in range(pos_integer):
        elements.append(i+1)
        elements.append((i+1) * -1)

    for perm in itertools.permutations(elements, pos_integer):
        if len(set([abs(x) for x in list(perm)])) == pos_integer:
            signed_perms.append(list(perm))

    return signed_perms

def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    pos_integer = parse_rosalind_input(args.input)

    signed_permutations = get_signed_permutations(pos_integer)

    print(len(signed_permutations))
    for signed_perm in signed_permutations:
        print(' '.join([str(perm) for perm in signed_perm]))

if __name__ == '__main__':
    main()