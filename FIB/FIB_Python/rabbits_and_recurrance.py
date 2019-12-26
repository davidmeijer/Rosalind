#!/usr/bin/env python3
""""

Author: David Meijer


Rosalind exercise: Rabbits and Recurrance

"""
import argparse

def define_arguments():
    """Parses input arguments.

    Function parses input arguments and returns an command line input
    object containing all input arguments.

    Args:
        command line arguments.

    Returns:
        parser (object): parsed command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')
    return parser

def parse_input(path_to_input):
    """ Parses input file.

    Parses assignment specific input file with two integers and returns
    both integers.

    Args:
        path_to_input (str): relative or absolute path to assignment
        specific input file.

    Returns:
        n (int): number of months after which a certain number of
        rabbits will be present.
        k (int): number of rabbit pairs in each litter.

    """
    with open(path_to_input, 'r') as in_fo:
        first_line = in_fo.readline().strip().split(' ')
        n, k = first_line[0], first_line[1]

    return int(n), int(k)

class Fibonacci():
    """For calculating sum of a Fibonacci recurrance relation.

    """
    def __init__(self, n, k, start=1):
        """

        Args:
            n (int): number of months after which a certain number of
            rabbit pairs will be present.
            k (int): number of rabbit pairs in each litter.
            maturity (int): number of months after which a generation
            becomes mature and starts reproducing.
            start (int): size of first generation
            fertility (int): number of pairs in litter.

        """
        self.n = n
        self.k = k
        self.start = start


    def pipe(self):
        """ Initiates full pipeline addressing all functions.

        """
        n2, n1 = self.set_start()

        for i in range(self.n):
            n2, n1 = self.next_generation(n2, n1)

        print(n1)

    def set_start(self):
        """Sets first two generations to initiate Fibonacci.

        Returns:
            n2 (int): number of rabbit pairs at start sequence (n-2).
            n1 (int): number of matured rabbit pairs at start sequence
            (n-1).

        """
        n2 = self.start
        n1 = self.start

        self.n = self.n - 2

        return n2, n1

    def next_generation(self, n2, n1):
        """Mature generation and generate offspring.

        Args:
            n2 (int): number of rabbit pairs at n-2.
            n1 (int): number of rabbit paris at n-1.

        Returns:
            new_n2 (int): new number of rabbit pairs at n-2.
            new_n1 (int): new number of rabbit pairs at n-1.

        """
        new_n2 = n1
        new_n1 = n1 + n2*self.k

        return new_n2, new_n1

def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    n, k = parse_input(args.input)
    Fibonacci(n, k).pipe()

if __name__ == '__main__':
    main()