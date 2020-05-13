#!/usr/bin/env python3
"""
Author: David Meijer

Script to solve Rosalind assignment Degree Array (deg).

Run script in command line as: "python3 deg.py <input_filename>"
"""

# import statements
from sys import argv

# functions
def parse_input(filename):
    """Parses edge list from input file in edge list format.

    filename: str, input filename to parse edge list from.

    Input file can only contain an edge list from a simple graph, one vertice
    per line.
    """
    edges = []

    with open(filename, 'r') as fo:
        for line in fo:
            edge = line.strip().split()
            edges.append(tuple(edge))

    return edges

def create_graph(edges):
    """Creates graph dictionary from edge list.

    edges: list of tuples, every tuple is an undirected vertice in a graph.
    graph: dict, keys as node and items as list of node neighbours.
    """
    graph = {}

    for edge in edges:
        node_a, node_b = edge[0], edge[1]

        if not node_a in graph:
            graph[node_a] = [node_b]
        else:
            graph[node_a].append(node_b)

        if not node_b in graph:
            graph[node_b] = [node_a]
        else:
            graph[node_b].append([node_a])

    return graph

# main
def main():

    # step 1: parse edge list from input file
    edges = parse_input(argv[1])

    # step 2: create graph from edge list
    graph = create_graph(edges)

    print(graph)

if __name__ == "__main__":
    main()
