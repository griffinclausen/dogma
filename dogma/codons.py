#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import product
from collections import defaultdict
import decimal
from random import choices

from dogma import (
    DEFAULT_DECIMAL_PRECISION,
    GeneticCode,
    DEFAULT_GENETIC_CODE,
    Nucleotide,
    DEGENERATE_NUCLEOTIDE_CODE,
    DEGENERATE_NUCLEOTIDE_CODE_REVERSED,
    is_valid_nucleotide_string,
    combine_nucleotides,
    STANDARD_CODONS,
    DEFAULT_CODON_LABEL
)

decimal.getcontext().prec = DEFAULT_DECIMAL_PRECISION


class Codon:
    """
    A sequence of three nucleotides.
    """

    def __init__(self, data=None, genetic_code=None):
        """
        Multiple input formats are acceptable
            n = Nucleotide('N')
            Codon([n,n,n]) == Codon((n, n, n))
        """

        if not isinstance(genetic_code, GeneticCode):
            genetic_code = DEFAULT_GENETIC_CODE
        self.genetic_code = genetic_code

        # aaa = Codon('AAA')
        # naa = Codon('NAA')
        # Xaa = combine_codons(aaa, naa)  # is dict
        # Codon(Xaa)
        if (isinstance(data, dict) and
                all([_ in STANDARD_CODONS for _ in data])):

            self.bases = None
            self.label = DEFAULT_CODON_LABEL
            self.members = list(data.keys())
            self.proportions = [data[_] for _ in self.members]

        # Codon(['AAA', 'GGG'])
        elif (isinstance(data, list) and
                all([isinstance(_, str) for _ in data]) and
                all([len(_) == 3 for _ in data]) and
                all([_ in STANDARD_CODONS for _ in data])):

            self.bases = None
            self.label = DEFAULT_CODON_LABEL
            self.members = list(set(data))
            self.proportions = [data.count(_) for _ in self.members]

        # Codon('AAA') or Codon('NNK') or ...
        # Codon([Nucleotide('A'),]*3) or ...
        # data as iterable of length 3 , or Nuc(),
        # or parameters that generates valid Nuc()
        elif (isinstance(data, (str, list, tuple)) and
                len(data) == 3):

            self.bases = [_
                          if isinstance(_, Nucleotide)
                          else Nucleotide(_)
                          for _ in data]
            self.label = ''.join([_.label for _ in self.bases])
            self.members = self.get_members()
            self.proportions = self.get_proportions()

        # Codon() --> Codon([Nucleotide(),]*3) or ...
        else:
            raise 'Error, Codon() parameters not valid'

        self.composition = dict(zip(self.members, self.proportions))

        self.degenerate = self.is_degenerate()
        self.equimolar = self.is_equimolar()

    def get_members(self):
        """
        Returns a list of codons (as triplet strings) with nonzero proportion.
        """
        b1, b2, b3 = [b.get_members() for b in self.bases]
        return [''.join(_) for _ in product(b1, b2, b3)]

    def get_proportions(self):
        """
        Returns a list of codon proportions, aligned with self.members
        """
        data = []
        for member in self.members:
            proportion = 1
            for m, b in zip(member, self.bases):
                proportion *= b.composition[m]
            data.append(proportion)
        return data

    def translate(self, dtype=None):
        """
        Returns a dictionary mapping amino acids (as single-letter codes)
        with proportions.
        """
        data = defaultdict(float)
        for c, p in zip(self.members, self.proportions):
            a = self.genetic_code[c]
            data[a] += p

        if dtype == 'decimal':
            data = {k: decimal.Decimal(v) for k, v in data.items()}
        elif dtype == 'float':
            data = {k: float(v) for k, v in data.items()}
        elif dtype == 'int':
            data = {k: int(v) for k, v in data.items()}
        elif dtype == 'object':
            pass
        elif dtype is None:
            data = {k: v for k, v in data.items()}
        else:
            # print('error', dtype, data)
            pass

        return data

    def samples(self, k=1):
        return [self.sample() for _ in range(k)]

    def sample(self):
        """
        Returns a single non-degenerate oligonucleotide string.
        """
        return choices(self.members, self.proportions)[0]

    def is_degenerate(self):
        """
        Codon is degenerate if there are multiple nonzero standard codons
        in composition.
        """
        return len(self.members) > 1

    def is_equimolar(self):
        """
        Codon is equimolar if all nonzero standard codons in composition have
        equal abundance.
        """
        return max(self.proportions) == min(self.proportions)

    def copy(self):
        return Codon([b.copy() for b in self.bases], self.genetic_code)

    def __str__(self):
        b0, b1, b2 = self.bases
        # return f'Codon(\n\t{b0}\n\t{b1}\n\t{b2}\n)'
        return f'Codon(\n\t{b0}\n\t{b1}\n\t{b2})'


def combine_codons(*codons, proportions=None):
    """
    codons unpacked to standard codons and proportions added

    ex: AAA + CCC + GGG + TTT --> [AAA, CCC, GGG, TTT]
    ex: NNN + AAA --> {AAA: 2 , AAC:1, AAG:1,.... TTT:1}
    """
    genetic_code = codons[0].genetic_code
    if proportions is None:
        proportions = [1, ] * len(codons)

    data = defaultdict(float)
    for codon, p in zip(codons, proportions):
        for sub_codon, sub_proportion in codon.composition.items():
            data[sub_codon] += p * sub_proportion

    return Codon(data, genetic_code)


def nucleotide_string_to_codons(s, output='Codon'):
    if is_valid_nucleotide_string(s):
        codon_strings = [s[3 * i:3 * i + 3] for i in range(len(s) // 3)]

        if output == 'Codon':
            return [Codon(*_) for _ in codon_strings]
        elif output == 'strings':
            return codon_strings


def degenerate_codon_string_to_standard_members(c):
    """
    Converts degenerate codon string to list of standard codons.
    """
    b1, b2, b3 = [DEGENERATE_NUCLEOTIDE_CODE[b] for b in c]
    return [''.join(_) for _ in product(b1, b2, b3)]


def combine_codon_labels(c1, c2):
    """
    Combines codon labels, assuming bases mixed randomly
    Ex: AAA + CCC + GGG + TTT -> NNN = (AAA, AAC, ....)
    Ex: ABC + CCC --> MBC
    """
    letters = []
    for b1, b2 in zip(c1, c2):
        b1 = b1.upper()
        b2 = b2.upper()
        _ = DEGENERATE_NUCLEOTIDE_CODE[b1] + DEGENERATE_NUCLEOTIDE_CODE[b2]
        base_letters = ''.join(sorted(list(set(_))))
        letters.append(DEGENERATE_NUCLEOTIDE_CODE_REVERSED[base_letters])
    return ''.join(letters)


def merge_codons(c1, c2, proportions=[1, 1], genetic_code=None):
    if genetic_code is None:
        genetic_code = DEFAULT_GENETIC_CODE
    nucs = [combine_nucleotides(n1, n2, proportions=proportions)
            for n1, n2 in zip(c1.bases, c2.bases)]
    return Codon(nucs, genetic_code)
