#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: Longest Increasing Subsequence
"""
import argparse
import copy

def define_arguments():
    """Defines possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind exercise specific input file.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parses Rosalind exericse specific input.

    Args:
        path_to_input (str): path to input file.

    Returns:
        n (int): positive integer.
        perm_n (list): list of integers containing a permutation of
        length n.

    """
    with open(path_to_input, 'r') as in_fo:
        n = in_fo.readline().strip()
        perm_n = in_fo.readline().strip().split()

    return int(n), [int(perm) for perm in perm_n]

def longest_subsequence(n, perm_list):
    """Get longest increasing subsequence of list of integers.

    """
    final_routes = []

    for j, seed in enumerate(perm_list):
        seed_list = []

        cut_perm_list = perm_list[j+1:]
        for i, next_seed in enumerate(cut_perm_list):
            if next_seed < seed:
                seed_list.append([seed, next_seed, cut_perm_list[i+1:]])

        pop_seeds = []
        for i, elem in enumerate(seed_list):
            if elem[-2] == 1 or len(elem[-1]) == 0 or elem[-2] < min(
                    elem[-1]):
                pop_seeds.append(i)
                final_routes.append(elem[:-1])
        for pop_seed in pop_seeds[::-1]:
            seed_list.pop(pop_seed)

        for seed_route in seed_list:
            for route in finalize_route(seed_route[:-1], seed_route[-1]):
                final_routes.append(route)

    return final_routes

def finalize_route(current_route, left_perms):
    """"""
    for left_perm in left_perms:
        if left_perm < current_route[-1]:
            new_route = current_route + [left_perm]
            next_perms = left_perms[(left_perms.index(left_perm))+1:]
            if new_route[-1] == 1 or len(next_perms) == 0 or new_route[-1] < min(next_perms):
                yield new_route
            else:
                for route in finalize_route(new_route, next_perms):
                    yield route

def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    n, perm_n = parse_rosalind_input(args.input)

    # Get longest increasing subsequence:
    incr_routes = longest_subsequence(n, perm_n[::-1])

    longest_routes = []
    longest_length = 0
    for route in incr_routes:
        if len(route) == longest_length:
            longest_routes.append(route)
        if len(route) > longest_length:
            longest_length = len(route)
            longest_routes = [route]

    #print(longest_routes)
    print(' '.join([str(x) for x in longest_routes[0][::-1]]))

    # Get longest decreasing subsequence:
    decr_routes = longest_subsequence(n, perm_n)

    longest_routes = []
    longest_length = 0
    for route in decr_routes:
        if len(route) == longest_length:
            longest_routes.append(route)
        if len(route) > longest_length:
            longest_length = len(route)
            longest_routes = [route]

    #print(longest_routes)
    print(' '.join([str(x) for x in longest_routes[0]]))

if __name__ == '__main__':
    main()