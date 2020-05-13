#!/usr/bin/env python3
"""
Author: David Meijer

Script for solving Rosalind assignment Connected Components (cc).

Run script in command-line as: "python3 cc.py <input_filename>"
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

def create_undirected_graph(n, edges):
    """Returns undirected graph as dictionary from node number and edge list.

    n: int, number of nodes in graph.
    edges: list of tuples as (a, b), a and b are nodes that are connected.
    graph: dict as {node : [neighbours], ...}, describes undirected graph.
    """
    graph = {}
    for node in range(1, n + 1):
        graph[node] = []

    for edge in edges:
        node_a, node_b = edge[0], edge[1]
        graph[node_a].append(node_b)
        graph[node_b].append(node_a)

    return graph

def connected_components(graph):
    """Returns number of connected components in graph.

    graph: dict as {node : [neighbours], ...}, describes undirected graph.
    cc: int, number of connected components.
    """
    unvisited = list(graph.keys())
    components = 0

    while unvisited:
        loc = unvisited[0]
        visited = explore(loc, graph)
        unvisited = list(set(unvisited) - set(visited))
        components += 1

    return components

def explore(root, graph):
    """Returns list of visitable nodes from root node.

    root: int, describing the starting node.
    graph: dict, describing directed graph as {node : [reachable nodes], ...}.
    visited: list, nodes that can visited from root.
    """
    visited = [root]
    queue = [root]

    while queue:
        loc = queue.pop(0)
        for neighbour in graph[loc]:
            if not neighbour in visited:
                queue.append(neighbour)
                visited.append(neighbour)

    return visited

# main
def main():

    # step 1: parse number of nodes n and edges m and edges from input file
    n, m, edges = parse_edge_list(argv[1])

    # step 2: create undirected graph from number of nodes and edge list
    graph = create_undirected_graph(n, edges)

    # step 3: calculate number of connected components in graph
    cc = connected_components(graph)

    print(cc)

if __name__ == "__main__":
    main()
