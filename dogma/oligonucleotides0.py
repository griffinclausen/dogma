#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dogma import (
    GeneticCode,
    DEFAULT_GENETIC_CODE,
    translate,
    Nucleotide,
    Codon,
)


class Oligonucleotide:
    """
    A sequence of nucleotides.
    """

    def __init__(self, data, genetic_code=None):
        """
        Multiple input formats are acceptable
            Oligonucleotide('AAANNK')
            Oligonucleotide(list('AAANNK'))
            Oligonucleotide(list(map(Nucleotide, list('AAANNK'))))
            Oligonucleotide([{'A':1}, {'A':0.5, 'C':0.5}]
            Oligonucleotide(list(map(Nucleotide, list('AAANNK'))))
        """

        if not isinstance(genetic_code, GeneticCode):
            genetic_code = DEFAULT_GENETIC_CODE
        self.genetic_code = genetic_code

        if isinstance(data, (str, list)):
            if isinstance(data[0], Codon):
                self.bases = []
                for c in data:
                    self.bases += c.bases
            else:
                self.bases = [_
                              if isinstance(_, Nucleotide)
                              else Nucleotide(_)
                              for _ in data]
        elif isinstance(data, Codon):
            self.bases = data.bases
            self.genetic_code = data.genetic_code
        else:
            self.bases = None

        self.length = self.get_length()
        self.label = self.get_label()
        self.compact_label = self.get_compact_label()

        self.codon_strings = self.get_codon_strings()
        self.codons = self.get_codons()

        self.get_base_profile()
        self.get_codon_profile()

    def get_bases(self):
        """
        Returns list of Nucleotides objects
        """
        return self.bases

    def get_codons(self):
        """
        Returns list of Codon objects derived from oligonucleotide bases
        """
        base_triplets = [self.bases[3 * i: 3 * i + 3]
                         for i in range(self.length // 3)]
        return [Codon(_, self.genetic_code) for _ in base_triplets]

    def get_length(self):
        """
        Returns integer length of oligonucleotide, number of bases
        """
        return len(self.bases)

    def get_label(self):
        """
        Returns a string of nucleotide labels.
        """
        return ''.join(_.label for _ in self.bases)

    def get_codon_strings(self):
        """
        Returns list of codon strings.
        """
        return [self.label[3 * i: 3 * i + 3] for i in range(self.length // 3)]

    def get_compact_label(self):
        """
        Returns shorthand label
        Attempts to reduce identical codon labels with numbers
        Ex: NNKNNKNNKAAANNK --> NNK<3>AAANNK
        """
        count = 1
        previous_string = ''
        output = []
        for c in self.get_codon_strings():
            if c == previous_string:
                count += 1
            else:
                if count > 1:
                    output.append(f'<{count}>')
                previous_string = c
                output.append(c)
                count = 1
        if count > 1:
            output.append(f'<{count}>')
        return ''.join(output)

    def translate(self, genetic_code=None):
        """
        Translates oligonucleotide into protein.

        Parameters:
        genetic_code: GeneticCode object
        """
        if not isinstance(genetic_code, GeneticCode):
            genetic_code = self.genetic_code
        return ''.join([genetic_code[_] for _ in self.get_codon_strings()])

    def samples(self, k=1, output='oligonucleotide', by_codon=True):
        return [self.sample(output=output, by_codon=by_codon)
                for _ in range(k)]

    def sample(self, output='oligonucleotide', by_codon=True):
        """
        Returns a single non-degenerate oligonucleotide string.

        If by_codon is True, sampling is based on the codons attribute,
        otherwise sampling is performed on each base individually.

        """
        if by_codon:
            oligonucleotide_string = ''.join(c.sample() for c in self.codons)
        else:
            oligonucleotide_string = ''.join([b.sample() for b in self.bases])

        if output == 'protein':
            return translate(oligonucleotide_string, self.genetic_code)
        elif output == 'both':
            return oligonucleotide_string, translate(oligonucleotide_string,
                                                     self.genetic_code)
        else:
            return oligonucleotide_string

    def get_base_profile(self):
        """
        Generates and returns base profile

        Usage
        -----
        >> nnk = Codon('NNK')
        >> oligo = Oligonucleotide(nnk)
        >> oligo.get_base_profile()
            {0: {'A':1, 'C':1 , 'G':1. 'T':1},
             1: {'A':1, 'C':1 , 'G':1. 'T':1},
             2: {               'G':1. 'T':1}}
        """
        self.base_profile = [b.composition for b in self.bases]
        return self.base_profile

    def __repr__(self):
        return f'Oligonucleotide({self.label})'


def reverse_complement(oligo):
    """
    """
    if isinstance(oligo, Oligonucleotide):
        oligo = oligo.label

    # TODO!


def get_nnk(n=1):
    """
    Simple helper function to easily generate NNK-based oligos.
    """
    supE = GeneticCode(1, {'TAG': 'Q'})
    oligo = Oligonucleotide('NNK' * n, supE)
    return oligo
