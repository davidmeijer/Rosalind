#!/usr/bin/env python3
"""

Author: David Meijer

Rosalind exercise: Completing a Tree.

"""
import argparse
import random

def define_arguments():
    """Defines possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')

    return parser

def parse_rosalinde_input(path_to_input):
    """Parses Rosalind exercise specific input.

    Args:
        path_to_input (str): path to input file.

    Returns:
        n (int): number of nodes in tree.
        edges (list): list as [[int, int], ...] containing known
        edges of tree.

    """
    edges = []
    with open(path_to_input, 'r') as in_fo:
        n = int(in_fo.readline().strip())
        for line in in_fo:
            edge = sorted(line.strip().split(' '))
            edges.append([int(x) for x in edge])

    return n, edges

def get_missing_nodes(n, edges):
    """Returns list of missing nodes.

    Args:
        n (int): number of nodes in tree.
        edges (list): list as [[int, int], ...] containing known edges
        of tree.

    Returns:
        nodes (list): list as [int, ...] containing all nodes
        missing from known edges.

    """
    # Nodes that should be present:
    nodes = [node+1 for node in range(n)]

    for edge in edges:
        for node in edge:
            if node in nodes:
                nodes.pop(nodes.index(node))

    return nodes

def connect_edges(edges):
    """Connects known edges into group.

    Args:
        edges (list): list as [[int, int], ...] containing known edges
        of tree.

    Returns:
        group (list): list as [[int, int], ...] containing known and
        grouped edges.
        edges (list): list as [[int, int], ...] containing left-over
        known edges after group creation.

    """
    group = []

    seed = random.choice(edges)
    group.append(seed)
    for member in group:
        for node in member:
            for edge in edges:
                if node in edge and not edge in group:
                    group.append(edge)

    for edge in group:
        edges.pop(edges.index(edge))

    return group, edges

def connect_branches(missing_nodes, branch_count):
    """Connects branches with missing nodes in least number of edges.

    Args:
        missing_nodes (list): list as [int, ...] containing all nodes
        missing from known edges.
        branch_count (int): number of groups with linked edges in tree.

    Returns:
        add_edges (int): minimum number of edges that can be added to
        the graph to produce a tree.

    """
    add_edges = branch_count + ( 1 * len(missing_nodes) ) - 1

    return add_edges

def main():
    """Main code.

    """
    # Parse Rosalind input:
    args = define_arguments().parse_args()
    n, edges = parse_rosalinde_input(args.input)

    # Get missing nodes:
    missing_nodes = get_missing_nodes(n, edges)

    # Get number of linked branches from tree:
    groups = []
    while len(edges) != 0:
        group, edges = connect_edges(edges)
        groups.append(group)
    branches = len(groups)

    # Connect branches with missing nodes in least number of edges:
    if len(missing_nodes) != 0:
        add_edges = connect_branches(missing_nodes, branches)
    else:
        add_edges = 0

    print(add_edges)

if __name__ == '__main__':
    main()