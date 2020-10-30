#!/usr/bin/env python3
"""
Author: David Meijer

Solution to Rosalind exercise 'Partial Sort'.
http://rosalind.info/problems/ps/
"""
from typing import List, Callable, Tuple, Optional
from argparse import ArgumentParser, Namespace
import numpy as np
from enum import Enum, unique
import sys


def define_arguments() -> Namespace:
    """
    Defines command line arguments and returns as argparse namespace.

    Out: args (Namespace) -- command line arguments.
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
    Out: k (int) -- k number of smallest elements to extract from arr.
    """
    with open(path, 'r') as fo:
        n = int(fo.readline().strip())
        arr = list(map(int, fo.readline().strip().split()))
        k = int(fo.readline().strip())
    return n, arr, k


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

        self.partial = None
        self.sorted = False

    def add(self, val: int):
        """
        Adds value to heap and heapifies to create extended heap.

        Arg: val (int) -- positive or negative integer to add to heap.
        """
        self._h[self._v] = val
        self._heapify(idx=self._v)
        self._v += 1

    def __str__(self) -> str:
        """
        Returns string representation of heap.

        Out: repr (str) -- returns heap as str.
        """
        if self.partial:
            return ' '.join(list(map(str, self._h[::-1][:self.partial])))

        return ' '.join(list(map(str, self._h)))

    def __repr__(self) -> Callable:
        """
        Returns method for representing heap.

        Out: __str__ (Callable) -- callable that returns string representation
            of heap.
        """
        return self.__str__()

    def sort(self, partial: Optional[int] = None) -> None:
        """
        Perform heap sort on heap.

        Arg: partial (int) -- stop sorting after getting `partial` steps
            (default: None).
        Out: sorted_heap (str) -- sorted heap.
        """
        self.partial = partial
        size = len(self._h)

        assert(partial <= size)

        if partial:
            stop_condition = size - partial
        else:
            stop_condition = 0

        while size > stop_condition:
            self._h[0], self._h[size - 1] = self._h[size - 1], self._h[0]
            size -= 1
            self._maintain_heap(idx=0, size=size)

        self.sorted = True

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

    def _maintain_heap(self, idx: int, size: int):
        """
        Maintains heap property by bubbling down node at `idx`.

        Arg: idx (int) -- node at position `idx` to bubble down.
        Arg: heap_size (int) -- current size unsorted part of heap.
        """
        l, r, parent = (2 * idx + 1), (2 * idx + 2), idx

        if self.type == self.HeapType.MAX:
            if l < size and self._h[l] > self._h[parent]:
                parent = l
            if r < size and self._h[r] > self._h[parent]:
                parent = r
            if parent != idx:
                self._h[idx], self._h[parent] = self._h[parent], self._h[idx]
                self._maintain_heap(idx=parent, size=size)

        elif self.type == self.HeapType.MIN:
            if l < size and self._h[l] < self._h[parent]:
                parent = l
            if r < size and self._h[r] < self._h[parent]:
                parent = r
            if parent != idx:
                self._h[idx], self._h[parent] = self._h[parent], self._h[idx]
                self._maintain_heap(idx=parent, size=size)


def create_heap(
        arr: List[int],
        heap_type: str,
        n: Optional[int] = None
    ) -> Heap:
    """
    Create a max binary heap from a list of integers.

    Arg: arr (list of int) -- list of positive or negative integers.
    Arg: n (int) -- length of arr (default: None).
    Arg: heap_type (str) -- type of binary heap to create
        (options: 'min' or 'max').
    Out: heap (Heap) -- Heap object containing max binary heap composed of
        values in arr.

    NOTE: if no n is given, n will be calculated by function from arr.
    """
    if not n:
        n = len(arr)
    heap = Heap(n, heap_type)
    for v in arr:
        heap.add(v)
    return heap


def main() -> None:
    """
    Driver code.
    """
    args = define_arguments()
    n, arr, k = parse_input(path=args.fi)
    heap = create_heap(arr=arr, n=n, heap_type='min')
    heap.sort(partial=k)
    print(heap)


if __name__ == "__main__":
    main()
