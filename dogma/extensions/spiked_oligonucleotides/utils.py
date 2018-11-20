from functools import reduce
import operator

from dogma import Nucleotide, rescale, STANDARD_NUCLEOTIDES


def prod(iterable):
    return reduce(operator.mul, iterable, 1)


def hamming_distance(a, b):
    """
    Counts number of nonequal matched elements in two iterables
    """
    return sum([int(i != j) for i, j in zip(a, b)])


def manhattan_distance(a, b):
    """
    Calculates the Manhattan (L1) distance between the values in two libraries.
    The union of the keys of both a and b are used, with the default value 0.
    """
    return sum([abs(a.get(k, 0) - b.get(k, 0)) for k in set().union(a, b)])


def combine_bases(b1, b2, proportions):
    """
    Adds two bases, b1 and b2, by adding weighted nucleotide compositions
    """
    data = {}
    for n in STANDARD_NUCLEOTIDES:
        data[n] = 0
        for nuc, p in zip([b1, b2], proportions):
            data[n] += rescale(nuc.composition).get(n, 0) * p
    return Nucleotide(rescale(data))
