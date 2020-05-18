#!/usr/bin/env python3
"""
Author:             David Meijer    
Description:        Scrip for solving Rosalind exercise Insertion Sort (ins).

Usage:              python3 ins.py <Rosalind input file>

Arguments:
    <Rosalind input file>   .txt input file containing a positive integer on
                            the first line and several space separated 
                            integers on the second line.
"""
# Import statements
from sys import argv

# Functions
def parse_input(filename):
    """Parses array of integers and the array length from the input file.
    
    Arguments:
        filename: str, input filename.
        
    Returns:
        array: list of ints, describing an unsorted array of integers.
        array_length: int, length of the array.
        
    The input file contains a positive integer on the first line and
    several sapce separated positive integers on the second line.
    """
    with open(filename, "r") as file_open:
        array_length = int(file_open.readline().strip())
        array = list(map(int, file_open.readline().strip().split()))
    
    return array_length, array
    
def insertion_sort(array):
    """Sorts an array through insert sort and returns the number of swaps.
    
    Arguments:
        array: list of ints, describing an unsorted array of integers.
        
    Returns:
        sorted_array: list of ints, describing a sorted array of integers.
        swaps: int, number of swaps performed for insert sorting the array.
    """
    swaps = 0

    for i in range(1, len(array)): 
  
        k = array[i] 
        j = i-1
        
        while j >= 0 and k < array[j] : 
                array[j+1] = array[j] 
                swaps += 1
                j -= 1
                
        array[j+1] = k
        
    return swaps, array
    
# Main
def main():
    
    # Define Rosalind input filename
    input_filename = argv[1]
    
    # Step 1: parse array and array length from Rosalind input file
    array_length, array = parse_input(input_filename)
    
    # Step 2: sort the array through insert sort and return the number of swaps
    swaps, sorted_array = insertion_sort2(array)
    
    # Step 3: print the number of swaps to stdout
    print(swaps)
    
if __name__ == "__main__":
    main()
    
