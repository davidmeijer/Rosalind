#!/usr/bin/env python3
"""

Author: David Meijer

Rosalind exercise: Calculating Expected Offspring.

"""
import argparse

def define_arguments():
    """Defines possible command line arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help='Rosalind input.')

    return parser

def parse_rosalind_input(path_to_input):
    """Parses Rosalind exericse specific input.

    Args:
        path_to_input (str): path to Rosalind exercise specific input
        file.

    Returns:
        AA_AA (int): number of couples with this genotype;
        AA_Aa (int): number of couples with this genotype;
        AA_aa (int): number of couples with this genotype;
        Aa_Aa (int): number of couples with this genotype;
        Aa_aa (int): number of couples with this genotype;
        aa_aa (int): number of couples with this genotype as a list in
        the same order.

    """
    with open(path_to_input, 'r') as in_fo:
        values = [int(x) for x in in_fo.readline().strip().split()]

    return values

def determine_dominant_phenotype(values, offspring=2):
    """Determine expected number of offspring displaying dom phenotype.

    Args:
        values (list): of integers --
        AA_AA (int): number of couples with this genotype;
        AA_Aa (int): number of couples with this genotype;
        AA_aa (int): number of couples with this genotype;
        Aa_Aa (int): number of couples with this genotype;
        Aa_aa (int): number of couples with this genotype;
        aa_aa (int): number of couples with this genotype as a list in
        the same order.
        offspring (int): number of offspring per couple.

    Returns:
        exp_offspring (float): expected number of offspring displaying
        dominant phenotype.

    """
    exp_offspring = []

    # AA-AA:
    pdom_AA_AA = 1
    exp_offspring.append(values[0] * pdom_AA_AA * offspring)

    # AA-Aa:
    pdom_AA_Aa = 1
    exp_offspring.append(values[1] * pdom_AA_Aa * offspring)

    # AA-aa:
    pdom_AA_aa = 1
    exp_offspring.append(values[2] * pdom_AA_aa * offspring)

    # Aa-Aa:
    pdom_Aa_Aa = 3/4
    exp_offspring.append(values[3] * pdom_Aa_Aa * offspring)

    # Aa-aa:
    pdom_Aa_aa = 1/2
    exp_offspring.append(values[4] * pdom_Aa_aa * offspring)

    # aa-aa:
    pdom_aa_aa = 0
    exp_offspring.append(values[5] * pdom_aa_aa * offspring)

    return sum(exp_offspring)

def main():
    """Main code.

    """
    args = define_arguments().parse_args()
    values = parse_rosalind_input(args.input)

    # Determine number of offspring with displaying dominant phenotype:
    exp_offspring = determine_dominant_phenotype(values)
    print(exp_offspring)

if __name__ == '__main__':
    main()