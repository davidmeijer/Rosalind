#!/usr/bin/env python3
"""
Author: David Meijer
Rosalind exercise: Introduction to Set Operations
"""
import argparse

def define_arguments():
    """Defines possible command line arguments.
    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', type=str, required=True,
                        help = 'Rosalind input file.')
                        
    return parser
    
def main():
    """Main code.
    
    """
    args = define_arguments().parse_args()
    
    
    with open(args.input, 'r') as fo:
        n = int(fo.readline().strip())
        A = set(list(map(int, fo.readline().strip()[1:-1].split(', '))))
        B = set(list(map(int, fo.readline().strip()[1:-1].split(', '))))
        
    U = set([i for i in range(1, n + 1)])
    
    with open('out_seto.txt', 'w') as fo:    
        fo.write(str(A | B) + '\n')
        fo.write(str(A & B) + '\n')
        fo.write(str(A - B) + '\n')
        fo.write(str(B - A) + '\n')
        fo.write(str(U - A) + '\n')
        fo.write(str(U - B) + '\n')

    
if __name__ == '__main__':
    main()
