#!/usr/bin/env python3
"""
Author: David Meijer

Solution to Dijkstra's Algorithm exercise on Rosalind.
25/10/2020
"""
from typing import Tuple, Dict, List, Union
from argparse import ArgumentParser, Namespace


WeightedGraph = Dict[int, Dict[int, int]]

def define_arguments() -> Namespace:
    """
    Returns command line arguments as Namespace object.

    Out: args (Namespace) -- command line aguments as Namespace.
    """
    parser = ArgumentParser()

    parser.add_argument(dest='input')

    parser.add_argument(
        '-output',
        dest='output',
        required=False,
        default='./dij_output.txt'
    )

    return parser.parse_args()


def parse_input(
    path: str
) -> Tuple[int, int, WeightedGraph]:
    """
    Returns Rosalind input edge list as hash table.

    Arg: path (str) -- path to input file.
    Out: n_n (int) -- number of nodes in graph.
    Out: n_e (int) -- number of edges in graph.
    Out: w_graph (dict) -- weighted graph as 
        {start (int) : {end (int) : weight (int), ...}, ...}.
    """
    with open(path, 'r') as fo:
        n_n, n_e = map(int, fo.readline().strip().split())
        w_graph = {n + 1 : {} for n in range(n_n)}
        for line in fo:
            start_n, end_n, w_e = map(int, line.strip().split())
            w_graph[start_n][end_n] = w_e

    return n_n, n_e, w_graph


def dijkstra(
    start_n: int,
    w_graph: WeightedGraph
) -> Dict[int, Union[int, float]]:
    """
    Returns lowest costs of start node to all other nodes in graph.

    Arg: start_n (int) -- start node in weighted graph.
    Arg: w_graph (dict) -- weighted graph as
        {graph (int) : {end (int) : weight (int), ...}, ...}.
    Out: costs (dict) -- dict as {end node (int) : cost, ...}.
    """
    costs = {
        n : (
            w_graph[start_n][n] 
            if n in w_graph[start_n] 
            else float('inf')
        ) for n in w_graph
    }

    costs[start_n] = 0

    parents = {
        n : (
            start_n
            if n in w_graph[start_n]
            else None
        )
        for n in w_graph
    }

    processed = []

    node = min(
        costs, 
        default=None, 
        key=costs.get
    )

    while node:
        cost = costs[node]
        nbs = w_graph[node]
        for nb in nbs:
            new_cost = cost + nbs[nb]
            if costs[nb] > new_cost:
                costs[nb] = new_cost
                parents[nb] = node
        processed.append(node)

        available = {
            node : cost for node, cost in costs.items()
            if node not in processed
        }

        node = min(
            available,
            default=None,
            key=available.get
        )

    return costs


def write_out_results(
    costs: Dict[int, Union[int, float]],
    path: str
) -> None:  
    """
    Writez out costs results to `path`.

    Arg: costs (dict) -- dict as {end node (int) : cost, ...}.
    Arg: path (str) -- output file path.

    NOTE: assumes dict keys are ordered (Python >=v3.7).
    """
    with open(path, 'w') as fo:
        fo.write(
            ' '.join([
                (
                    str(item) if not item == float('inf') else str(-1)
                ) 
                for k, item in costs.items()
            ])
        )


def main() -> None:
    """
    Driver code.
    """
    args = define_arguments()

    _, _, w_graph = parse_input(path=args.input)

    costs = dijkstra(
        start_n=1,
        w_graph=w_graph
    )

    write_out_results(
        costs=costs,
        path=args.output
    )

if __name__ == '__main__':
    main()