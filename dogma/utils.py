#!/usr/bin/env python
# -*-encoding=utf-8-_*-


"""Package-wide variables and helper functions."""

from itertools import product
from random import choice


DEFAULT_DECIMAL_PRECISION = 200

STANDARD_NUCLEOTIDES = list('ACGT')

DEGENERATE_NUCLEOTIDE_CODE = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': 'AG',
    'Y': 'CT',
    'S': 'CG',
    'W': 'AT',
    'K': 'GT',
    'M': 'AC',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG',
    'N': 'ACGT'
}

DEGENERATE_NUCLEOTIDES = sorted(DEGENERATE_NUCLEOTIDE_CODE.keys())

DEGENERATE_NUCLEOTIDE_CODE_REVERSED = {v: k for k, v in
                                       DEGENERATE_NUCLEOTIDE_CODE.items()}

DEGENERATE_NUCLEOTIDE_CODE_COMPOSITION = {code: {n: int(n in _)
                                          for n in STANDARD_NUCLEOTIDES}
                                          for code, _ in
                                          DEGENERATE_NUCLEOTIDE_CODE.items()}

NUCLEOTIDE_BASE_PAIRS = dict(zip('ACGT', 'TGCA'))

_D = DEGENERATE_NUCLEOTIDE_CODE_REVERSED
DEGENERATE_NUCLEOTIDE_PAIRS = {k: _D[''.join(sorted(
                                     [NUCLEOTIDE_BASE_PAIRS[_]
                                      for _ in v]))]
                               for k, v in DEGENERATE_NUCLEOTIDE_CODE.items()}

DEFAULT_NUCLEOTIDE_LABEL = 'N'

DEGENERATE_NUCLEOTIDE_PAIRS_ = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'R': 'Y',
    'Y': 'R',
    'S': 'S',
    'W': 'W',
    'K': 'M',
    'M': 'K',
    'B': 'V',
    'D': 'H',
    'H': 'D',
    'V': 'B',
    'N': 'N'
}


DEFAULT_OLIGONUCLEOTIDE_LABEL = 'NNN'

STANDARD_CODONS = [''.join(_) for _ in
                   product(STANDARD_NUCLEOTIDES, repeat=3)]

DEGENERATE_CODONS = [''.join(_) for _ in
                     product(DEGENERATE_NUCLEOTIDES, repeat=3)]

DEFAULT_CODON_LABEL = 'NNN'

STANDARD_AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

DEFAULT_AMINO_ACID_LABEL = 'X'

STOP_LABEL = '*'

DEFAULT_RESIDUE_LABEL = 'X'


def rescale(data, total=1):
    """
    Rescales numerical values in lists or dictionary values to sum
    to specified total.

    Usage
    *****
    rescale([1, 3]) -> [0.25 0.75]
    rescale({'a': 1, 'b':'9']) -> {'a': 0.1, 'b': 0.9}
    """

    if isinstance(data, list):
        input_total = sum(data)
        assert input_total != 0, 'Error in dogma.rescale(), input_total == 0'

        return [_ / input_total * total for _ in data]

    elif isinstance(data, dict):
        input_total = sum(data.values())
        assert input_total != 0, 'Error in doe.rescale(), input_total == 0'

        return {k: v / input_total * total for k, v in data.items()}
    else:
        return None


def get_frequency_dictionary(data):
    """
    Takes a string or list of strings and returns a dictionary of unique
    members and their abundance. Similar to itertools.Counter.
    """
    return {k: data.count(k) for k in set(data)}


def get_random_oligonucleotide(length=3, letters='ACGTRYSWKMBDHVN'):
    """
    Returns string of random degenerate nucleotides.

    Parameters
    **********
    length - int - length of string to return
    letters - string or list - collection of letters to select from, with replacement.

    Each of the 15 letters 'ACGTRYSWKMBDHVN' are equally likely.
    Alternatively, use letters='ACGT' for standard oligo strings
    """
    return ''.join(choice(list(letters)) for _ in range(length))
