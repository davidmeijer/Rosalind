#!/usr/bin/env python3
"""

Author: David Meijer

Rosalind exercise: Enumerating k-mers Lexicographically.

"""
import argparse
import itertools

def define_arguments():
    """Defines possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parse exercise specific Rosalind input.

    Args:
        path_to_input (str): relative or absolute path to exercise
        specific Rosalind input file.

    Returns:
        alphabet (list): containing ordered symbols.
        n (int): positive integer defining string length for strings
        that can be formed from the alphabet.

    """
    with open(path_to_input, 'r') as in_fo:
        alphabet =  in_fo.readline().strip().split(' ')
        n = int(in_fo.readline().strip())

    return alphabet, n

def get_permutations(elem_list, n):
    """Return all permutations of element list.

    Args:
        elem_list (list): list of lexicographically ordered items.
        n (int): positive integer defining string length for strings
        that can be formed from the alphabet.

    Returns:

    """
    elem_list_permutations = []

    for perm_tup in itertools.permutations(elem_list, n):
        elem_list_permutations.append(''.join(perm_tup))

    return list(set(elem_list_permutations))

def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    alphabet, n = parse_rosalind_input(args.input)

    # Multiply the alphabet by n to also include permutations later on
    # with itself:
    new_alphabet =  alphabet*n

    total_permutations = []
    for i in range(n):
        alphabet_permutations = get_permutations(new_alphabet, i+1)
        for perm in alphabet_permutations:
            total_permutations.append(perm)

    # Create new order:
    new_order = ''.join(alphabet)

    # Sort kmers with new lexicographical order:
    ordered_permutations = sorted(total_permutations,
                    key=lambda word: [new_order.index(c) for c in word])

    for kmer in ordered_permutations:
        print(kmer)

if __name__ == '__main__':
    main()