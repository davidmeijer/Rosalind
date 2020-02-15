#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: Counting Subsets
"""
import argparse

def define_arguments():
    """Defines possible command line arguments.
    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help = 'Rosalind exercise input file.')
                        
    return parser

def nCr(n, r):
    """Computes number of unique picks from a set.
    
    Args:
        n (int): total size of set.
        r (int): size of pick.
    
    Returns:
        k (int): number of unique picks of size r from n.
    
    """
    if r > 0 and r != n:
        k = factorial(n) / (factorial(n - r) * factorial(r))
    else: k = 1

    return int(k) % 1000000    

def factorial(x):
    """Computes factorial (!) of x.
    
    Args:
        x (int): integer for which factorial is computed.
        
    Returns:
        x (int): factorial of x.
    
    """
    
    for j in [i for i in range(x - 1, 0, -1)]:
        x *= j
    
    return x
    
def main():
    """Main code.
    
    """
    args = define_arguments().parse_args()
    
    with open(args.input, 'r') as fo: n = int(fo.readline().strip())
    
    nsubset = 0
    for r in [i for i in range(n + 1)]: nsubset += nCr(n, r)
    
    answer1 = nsubset % 1000000 (nCr)
    
    answer2 = (2 ** n) % 1000000 (2 ** r)
    # You can choose ON/OFF (2 things), and we choose r of them.
    
    print(answer1, answer2, n)
        
    
        
if __name__ == '__main__':
    main()
