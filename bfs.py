#!/usr/bin/env python3
"""
Author: David Meijer

Script for solving Rosalind assignment Breadth-First Search (bfs).

Run script in command-line as: "python3 bfs.py <input_filename>"
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

def simple_directed_graph(n, edges):
    """Returns a simple directed graph made from edge list as dictionary.

    n: int, number of nodes.
    edges: list of tuples as (a, b), describing edge between nodes a and b.
    """
    nodes = [node for node in range(1, n + 1)]

    graph = {}
    for node in nodes:
        graph[node] = []

    for edge in edges:
        start, destination = edge[0], edge[1]
        graph[start].append(destination)

    return graph

def shortest_paths(root, graph):
    """Returns list of distances of root towards all other nodes in graph.

    root: int, describing the starting node.
    graph: dict, describing directed graph as {node : [reachable nodes], ...}.
    """
    distances = {}
    for node in graph.keys():
        distances[node] = -1

    distances[root] = 0
    queue = [root]

    while queue:
        loc = queue.pop(0)
        for neighbour in graph[loc]:
            if distances[neighbour] == -1:
                queue.append(neighbour)
                distances[neighbour] = distances[loc] + 1

    return distances

# main
def main():

    # step 1: parse number of nodes n and edges m and edges from input file
    n, m, edges = parse_edge_list(argv[1])
    # m remains unused

    # step 2: create simple directed graph from nodes n and edges
    graph = simple_directed_graph(n, edges)

    # step 3: calculates lengths of shortest paths between node 1 and all other
    distances = shortest_paths(1, graph)

    # step 4: print results to stdout
    int_result = [distances[node] for node in sorted(graph.keys())]
    str_result = list(map(str, int_result))
    print(" ".join(str_result))

if __name__ == "__main__":
    main()
