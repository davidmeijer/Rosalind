#!/usr/bin/env python3
"""

Author: David Meijer

Rosalind exercise: Genome Assembly as Shortest Superstring.

"""
import argparse
import copy

def define_arguments():
    """Defines possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parses exercise specifc Rosalind input.

    Args:
        path_to_input (str): path to exercise specific input file.

    Returns:
        fasta_dict (dict): dictionary containing parses fasta input
        file as {header : seq, ...}.

    """
    fasta_dict = {}

    with open(path_to_input, 'r') as in_fo:
        for line in in_fo:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                fasta_dict[header] = []
                continue
            else:
                fasta_dict[header].append(line)

    for header, seq in fasta_dict.items():
        fasta_dict[header] = ''.join(seq)

    return fasta_dict

def partition_ends(seq):
    """Partition end parts of sequence untill half.

    Args:
        seq (str): DNA, RNA or peptide sequence.

    Returns:
        start_snippets (list): list of strings containing start snippets
        for genome assembly.
        end_snippets (list): list of strings containing end snippets
        for genome assembly.

    """
    start_snippets = []
    end_snippets = []

    for pos in range(int(0.5 * len(seq))):
        pos += int(0.5 * len(seq))
        # Get snippets first half:
        start_snippets.append(seq[:pos+1])
        # Get snippets end half:
        end_pos = len(seq) - pos - 1
        end_snippets.append(seq[end_pos:])

    return start_snippets, end_snippets

def match_snippets(snippets_dict):
    """Matches start with end snippets and vica versa.

    Args:
        snippets_dict (dict): contains header with start and end
        snippets.

    Returns:
        matching_dict (dict): contains header with directional
        matches as {header : {start: header_match, end: header_match},
        ...}.

    """
    matching_dict = {}

    for header in snippets_dict:
        matching_dict[header] = {}
        # Match starting snippets with ending snippets:
        # Make sure longest matches are stored first in the matching_dict:
        matching_dict[header]['start'] = []
        for snippet in snippets_dict[header]['start'][::-1]:
            for comp_header in snippets_dict:
                if snippet in snippets_dict[comp_header]['end']:
                    if header != comp_header:
                        matching_dict[header]['start'].append([comp_header, snippet])

        matching_dict[header]['end'] = []
        for snippet in snippets_dict[header]['end'][::-1]:
            for comp_header in snippets_dict:
                if snippet in snippets_dict[comp_header]['start']:
                    if header != comp_header:
                        matching_dict[header]['end'].append([comp_header, snippet])

    return matching_dict

class pathfinding():
    """Finds possible assembly paths.

    """
    def __init__(self, matching_dict):
        """Assigns variables.

        Args:
            matching_dict (dict): contains header with directional
            matches as {header : {start: header_match, end: header_match},
            ...}.

        """
        self.matching_dict = matching_dict
        self.routes = []

    def find_routes(self):
        """Filters out short matches when there is a longer match.

        """
        # Maybe first determine from which seq you can go through all other
        # sequences? Best final sequence has most overlap between sequences
        # resulting in the shortest assembled genome!

        # Restructure below so it takes every sequence and index every
        # possible route per starting position. Than determine starting
        # position.

        for header, matches in list(self.matching_dict.items()):
            route_dict = copy.deepcopy(self.matching_dict)
            del route_dict[header]
            paths = self.find_paths(header, matches, route_dict)
            for path in paths:
                self.routes.append(path)

        return self.routes

    def find_paths(self, header, matches, route_dict):
        """Find assembly path through sequences.

        Args:
            header (str): header of start sequence.
            matches (dict): overlapping sequences at end and beginning
            of sequence as {'start' : [[overlapping_header, match], ...],
            'end' : [[overlapping_header, match], ...]}
            route_dict (dict): contains header with directional
            matches as {header : {start: header_match, end: header_match},
            ...}. Start sequence is deleted from this dictionary.

        Returns:
            paths ():

        """
        paths = []

        # Add header and first steps as two-length paths:
        for match in matches['end']:
            current_route = [header]
            current_route.append(match[0])
            paths.append(current_route)

        path_lengths = [0]
        iterations = 0
        while iterations < len(list(self.matching_dict.keys())):
            paths = self.next_step(route_dict, paths)
            for path in paths:
                path_lengths.append(len(path))
            iterations += 1

        return paths

    def next_step(self, route_dict, paths):
        """Take next step in pathfinding.

        """
        new_paths = []

        for path in paths:
            last_seq = path[-1]
            next_steps = route_dict[last_seq]['end']
            for next_step in next_steps:
                path = copy.deepcopy(path)
                if not next_step[0] in path:
                    path.append(next_step[0])
            new_paths.append(path)

        return new_paths

def filter_routes_on_length(routes, filter_length):
    """Filters out routes that do not contain certain seq length.

    Args:
        routes (list): list of sequences headers for certain route.
        filter_length (int): length that needs to be filtered on.

    Returns:
        routes (list): all routes with certain filter_length.

    """
    new_routes = []

    for route in routes:
        if len(route) == filter_length:
            new_routes.append(route)

    return new_routes

def assemble_route(route, fasta_dict, matching_dict):
    """Assemble route with snippets and fasta dicts.

    Args:
        route (list):
        fasta_dict (dict):
        matching_dict (dict):

    Returns:
        assembly (str):

    """
    assembly = []

    # Start assembly:
    assembly.append(fasta_dict[route[0]])

    # Elongate assembly:
    for i, header in enumerate(route[1:]):
        seq = fasta_dict[header]

        for match in matching_dict[route[i]]['end']:
            if match[0] == header:
                overlap = match[1]

        non_overlapping_seq = seq.partition(overlap)[2]
        assembly.append(non_overlapping_seq)

    return ''.join(assembly)

def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    fasta_dict = parse_rosalind_input(args.input)

    snippets_dict = {}
    for header, seq in fasta_dict.items():
        snippets_dict[header] = {}
        start_snippets, end_snippets = partition_ends(seq)
        snippets_dict[header]['start'] = start_snippets
        snippets_dict[header]['end'] = end_snippets

    # Start snippets can only be matched with end snippets and vica versa:
    matching_dict = match_snippets(snippets_dict)
    routes = pathfinding(matching_dict).find_routes()
    routes = filter_routes_on_length(routes, len(list(fasta_dict.keys())))

    # Assemble every route:
    shortest_assembly_length = float('inf')
    shortest_assembly = []
    for route in routes:
        assembly = assemble_route(route, fasta_dict, matching_dict)
        if len(assembly) < shortest_assembly_length:
            shortest_assembly = [assembly]
            shortest_assembly_length = len(assembly)
        elif len(assembly) == shortest_assembly_length:
            shortest_assembly.append(assembly)

    print(shortest_assembly[0])

if __name__ == '__main__':
    main()