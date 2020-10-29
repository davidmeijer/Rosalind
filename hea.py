#!/usr/bin/env python3
"""
Author: David Meijer

Solution to Rosalind exercise 'Building a Heap'.
http://rosalind.info/problems/hea/
"""
from typing import Tuple, List, Callable, Optional
from argparse import ArgumentParser, Namespace
import numpy as np
from enum import Enum, unique
import sys


def define_arguments() -> Namespace:
    """
    Defines possible command line arguments.

    Out: args (Namespace) -- user defined CL arguments.
    """
    parser = ArgumentParser()
    parser.add_argument(dest='fi')
    return parser.parse_args()


def parse_input(path: str) -> Tuple[int, List[int]]:
    """
    Parses input data from Rosalind input file.

    Arg: path (str) -- path to Rosalind input file.
    Out: n (int) -- number of nodes.
    Out: arr (list of int) -- node values.
    """
    with open(path, 'r') as fo:
        n = int(fo.readline().strip())
        arr = list(map(int, fo.readline().strip().split()))
    return n, arr


class Heap:
    """
    Stores heap types and heap.
    """
    @unique
    class HeapType(Enum):
        """
        Stores binary max and binary min heap types.
        """
        MAX = 1
        MIN = 2

    def __init__(self, n: int, heap_type: str) -> None:
        """
        Initializes the binary max or binary min heap.

        Arg: n (int) -- length of the heap.
        Arg: heap_type (str) -- type of binary heap to create
            (options: 'min' or 'max').
        """
        self._h = np.empty(n, dtype=object)
        self._v = 0

        try:
            self.type = self.HeapType[heap_type.upper()]
        except KeyError:
            print(f'Invalid heap type: {heap_type}')
            sys.exit(1)

    def add(self, val: int):
        """
        Adds value to heap and heapifies to create extended heap.

        Arg: val (int) -- positive or negative integer to add to heap.
        """
        self._h[self._v] = val
        self._heapify(idx=self._v)
        self._v += 1

    def _heapify(self, idx: int):
        """
        Heapify the heap with new value at index `idx`.

        Arg: idx (int) -- positive integer pointer to newly added value to
            heap.
        """
        if idx == 0:
            return

        parent = (idx - 1) // 2
        child = idx

        if (self.type == self.HeapType.MAX and
                self._h[parent] > self._h[child]):
            return
        elif (self.type == self.HeapType.MIN and
                self._h[parent] < self._h[child]):
            return
        else:
            self._h[parent], self._h[child] = self._h[child], self._h[parent]
            self._heapify(parent)

    def __str__(self) -> str:
        """
        Returns string representation of heap.

        Out: repr (str) -- returns heap as str.
        """
        return ' '.join(list(map(str, self._h)))

    def __repr__(self) -> Callable:
        """
        Returns method for representing heap.

        Out: __str__ (Callable) -- callable that returns string representation
            of heap.
        """
        return self.__str__()


def create_max_heap(arr: List[int], n: Optional[int] = None) -> Heap:
    """
    Create a max binary heap from a list of integers.

    Arg: arr (list of int) -- list of positive or negative integers.
    Arg: n (int) -- length of arr (default: None).
    Out: heap (Heap) -- Heap object containing max binary heap composed of
        values in arr.

    NOTE: if no n is given, n will be calculated by function from arr.
    """
    if not n:
        n = len(arr)
    heap = Heap(n, 'max')
    for v in arr:
        heap.add(v)
    return heap


def main() -> None:
    """
    Driver code.
    """
    args = define_arguments()
    n, arr = parse_input(path=args.fi)
    heap = create_max_heap(n=n, arr=arr)
    print(heap)


if __name__ == '__main__':
    main()
