#!/usr/bin/env python3
"""
Author: David Meijer

Script for solving Rosalind assignment Double-Degree Array (ddeg).

Run script in command-line as: "python3 ddeg.py <input_filename>"
"""

# import statements
from sys import argv

# functions
def parse_edge_list(filename):
    """Parses number of nodes and edges and edges from input edge list file.

    filename: str, input file in edge list format.

    File in edge list format has number od nodes and edges on the first line
    and the edges on consecutive lines.
    """
    edges = []
    with open(filename, 'r') as fo:
        n, m = fo.readline().strip().split()
        for line in fo:
            edge = map(int, line.strip().split())
            edges.append(tuple(edge))

    return int(n), int(m), edges

def simple_graph(n, m, edges):
    """Returns a simple graph made from edge list as dictionary.

    n: int, number of nodes.
    m: int, number of edges.
    edges: list of tuples as (a, b), describing edge between nodes a and b.
    """
    nodes = [node for node in range(1, n + 1)]

    graph = {}
    for node in nodes:
        graph[node] = []

    for edge in edges:
        node_a, node_b = edge[0], edge[1]
        graph[node_a].append(node_b)
        graph[node_b].append(node_a)

    return graph

def neighbour_degrees(node, graph):
    """Returns the sum of degrees of all node's neighbours.

    node: int, describes node in graph.
    graph: dict as {node : [neighbour node, ...]}, describing graph structure.
    degree: int, sum of degrees of all node's neighbours.
    """
    degrees = []
    for neighbour in graph[node]:
        degrees.append(len(graph[neighbour]))

    return sum(degrees)

# main
def main():

    # step 1: parse number of nodes n and edges m and edges from input file
    n, m, edges = parse_edge_list(argv[1])

    # step 2: create undirected graph from
    graph = simple_graph(n, m, edges)

    # step 3: calculate the degree of every node
    degrees = []
    for node in sorted(graph.keys()):
        degrees.append(neighbour_degrees(node, graph))

    print(" ".join(list(map(str, degrees))))

if __name__ == "__main__":
    main()
