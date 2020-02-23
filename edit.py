#!/usr/bin/env python3
"""Author: David Meijer"""

import argparse

def DefineArguments():
    """Defines possible command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type = str, required = True,
                        help = 'Rosalind exercise input file.')

    return parser

def FastaParser(fn):
    """Parser FASTA file into dictionary."""
    fasta = {}

    with open(fn, 'r') as fo:
        for line in [line.strip() for line in fo]:
            if line.startswith('>'):
                header = line
                fasta[header] = ''
            else:
                fasta[header] += line

    return fasta

class Matrix(object):
    """Matrix class."""
    def __init__(self, nrows = 0, ncols = 0, fill = 0):
        self._nrows = nrows
        self._ncols = ncols
        self._fill = fill

        self.SetMatrix()

    def __str__(self):
        rows = []

        for row in self.matrix:
            rows.append(str(' '.join([str(item) for item in row])))

        return '\n'.join(rows)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.matrix)

    @property
    def nrows(self):
        return self._nrows

    @property
    def ncols(self):
        return self._ncols

    @property
    def fill(self):
        return self._fill

    @nrows.setter
    def nrows(self, new_nrows):
        self._nrows = new_nrows
        self.SetMatrix()

    @ncols.setter
    def ncols(self, new_ncols):
        self._ncols = new_ncols
        self.SetMatrix()

    @ncols.setter
    def fill(self, new_fill):
        self._fill = new_fill
        self.SetMatrix()

    def SetMatrix(self):
        """Sets matrix with defined values."""
        self.matrix = [[self._fill for i in range(self.ncols)] \
                        for j in range(self.nrows)]

    def SetValue(self, x, y, value):
        """Sets value in matrix at intersection of row x and column y."""
        self.matrix[x][y] = value
        return self.matrix

    def ReturnValue(self, x, y):
        """Returns value in matrix at intersection of row x and column y."""
        return self.matrix[x][y]

    def Dimensions(self):
        """Return matrix dimensions."""
        return [self.ncols, self.nrows]

    def SetDiagonal(self, value):
        """Sets diagonal of matrix to certain value."""

    def Transpose(self):
        """Return transposed matrix."""

    def MultiplyMatrix(self, other):
        """Returns product of matrices."""

    def SubtractMatrix(self, other):
        """Returns subtraction of matrices."""

    def AddMatrix(self, other):
        """Returns addition of matrices."""


class Protein(object):
    """Protein class."""
    def __init__(self, string):
        self.string = string

    def __str__(self):
        return self.string

    def __repr__(self):
        return self.string

    def __eq__(self, other):
        return self.string == other.string

    def __hash__(self):
        return hash(self.string)

    def Length(self):
        """Returns number of amino acids in peptide string."""
        return len(self.string)

    def LevenshteinDistance(self, other):
        """Calculates Levenshtein distance between peptides."""
        dist = Matrix(self.Length() + 1, other.Length() + 1)

        for i in range(0, self.Length() + 1):
            dist.SetValue(i, 0, i)

        for j in range(0, other.Length() + 1):
            dist.SetValue(0, j, j)

        for j in range(1, other.Length() + 1):
            for i in range(1, self.Length() + 1):
                if self.string[i - 1] == other.string[j - 1]:
                    substitutionCost = 0
                else:
                    substitutionCost = 1

                dist.SetValue(i, j, min([
                dist.ReturnValue(i - 1, j) + 1,
                dist.ReturnValue(i, j - 1) + 1,
                dist.ReturnValue(i - 1, j - 1) + substitutionCost
                ]))

        return dist.ReturnValue(self.Length(), other.Length())

    def CommonLongestSubstring(self, other):
        """Calculates common longest substring between peptides."""


def main():
    args = DefineArguments().parse_args()
    fasta = FastaParser(args.i)

    t = Protein(fasta[list(fasta.keys())[0]]) # Ugh...
    s = Protein(fasta[list(fasta.keys())[1]]) # Ugh...

    print(s.LevenshteinDistance(t))

if __name__ == '__main__':
    main()
