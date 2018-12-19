#!/usr/bin/env python
# ~*~# -*- coding: utf-8 -*-

"""
The genetic code defines the mapping between nucleic acids and amino acids.

This module defines the GeneticCode object and related functions.

Usage
*****
std = GeneticCode()  # standard genetic code
std['AAA'] -> 'K'
std['BAD_CODON'] -> '*' # stop symbol
"""


from dogma import (
    NBCI_GENETIC_CODE_NAMES,
    NCBI_CODONS,
    NCBI_AMINO_ACIDS
)


class GeneticCode:
    """
    System for translating codons into amino acids.

    Parameters
    ----------
    ncbi_id: integer witin NBCI_GENETIC_CODE_NAMES keys,
             default is 1, standard genetic code
    updated_mappings: dictionary
    stop_symbol: string with length 1
    error_symbol: string with length 1
    name: str
    """

    def __init__(self,
                 ncbi_id=1,
                 updated_mappings={},
                 stop_symbol='*',
                 error_symbol='_',
                 name=''):

        if ncbi_id in NBCI_GENETIC_CODE_NAMES:
            name = NBCI_GENETIC_CODE_NAMES[ncbi_id]
            codons = NCBI_CODONS
            amino_acids = NCBI_AMINO_ACIDS[ncbi_id]
        else:
            codons = []
            amino_acids = ''

        if updated_mappings:
            assert isinstance(updated_mappings, dict)
            assert all([(len(c) == 3 and isinstance(c, str))
                        for c in updated_mappings.keys()])
            assert all([(len(a) == 1 and isinstance(a, str))
                        for a in updated_mappings.values()])
            updated_mappings = {k.upper().replace('U', 'T'): v
                                for k, v in updated_mappings.items()}

        if error_symbol is not '_':
            assert isinstance(error_symbol, str) and len(error_symbol) == 1
        self.error_symbol = error_symbol

        if stop_symbol is not '*':
            assert isinstance(stop_symbol, str) and len(stop_symbol) == 1
            amino_acids = amino_acids.replace('*', stop_symbol)
            self.stop_symbol = stop_symbol

        if isinstance(name, str):
            self.name = name

        # define genetic code mapping from codons to amino_acids
        self.code = dict(zip(codons, amino_acids))

        # overwrite specified codon to amino acid mappings
        for codon, amino_acid in updated_mappings.items():
            self.code[codon] = amino_acid

        # create list of codons and amino acids
        self.codons = codons
        self.amino_acids = amino_acids

        # define reversed genetic code, where each amino acids
        # maps to a list of codons
        self.code_reversed = {a: [_ for _ in codons
                                  if self.code[_] == a]
                              for a in amino_acids}

    def update(self, updated_mappings={}, expand=False):
        """
        Updates genetic code, reprogramming codons to amino acids

        Parameters
        ----------
        updated_mappings: dict

        expand: Boolean
            if True, genetic code is 'expanded' to include new codons
            if False, only keys in GeneticCode.code.keys() will be remapped,
                      others silently ignored.
            Default is False

        Usage
        -----
        standard = GeneticCode()  #  default parameter ncbi_id = 1
                                  #  corresponds to 'Standard' code
        supE = GeneticCode(updated_mappings={'TAG': 'Q'})
        """
        assert isinstance(updated_mappings, dict)
        assert all([(len(c) == 3 and isinstance(c, str))
                    for c in updated_mappings.keys()])
        assert all([(len(a) == 1 and isinstance(a, str))
                    for a in updated_mappings.values()])

        if not expand:
            updated_mappings = {k: v
                                for k, v in updated_mappings.items()
                                if k in self.code}

        for codon, amino_acid in updated_mappings.items():
            self.code[codon] = amino_acid

    def get_codons(self):
        """
        Returns sorted list of codons in genetic code.
        """
        return sorted(list(self.code.keys()))

    def get_amino_acids(self):
        """
        Returns list of amino acids, ordered by sorted codons.
        """
        return [self.code[codon] for codon in self.get_codons()]

    def __getitem__(self, codon):
        """
        Returns value of code dictionary for specified codon key.

        Parameter
        ---------
        codon: string within code dictionary

        Returns
        -------
        One-letter amino acid string

        Usage
        -----
        standard = GeneticCode(1)
        standard['GGG'] --> 'G'
        """
        codon = codon.upper().replace('U', 'T')

        return self.code.get(codon, self.error_symbol)

    def __setitem__(self, codon, amino_acid):
        """
        Sets value of code dictionary for specified codon key
        to specified amino acid.

        Parameter
        ---------
        codon: string within genetic code dictionary

        Usage
        -----
        >> std = GeneticCode()
        >> std['GGG'] --> 'G'
        >> std['GGG'] = 'Z'
        >> std['GGG'] --> 'Z'

        """
        codon = codon.upper().replace('U', 'T')
        if codon not in self.code:
            raise "Attempting to set codon not in genetic code."
        else:
            self.code[codon] = amino_acid

    def __str__(self):
        """
        Returns readable genetic code list
        """
        codons = self.get_codons()
        amino_acids = self.get_amino_acids()
        return '\n'.join([' '.join([c, a])
                         for c, a in zip(codons, amino_acids)])

    def __repr__(self):
        """
        Returns string that would evaluate to an equivalent
        GeneticCode instance.
        """
        codons = self.get_codons()
        amino_acids = self.get_amino_acids()
        return f'GeneticCode(codons={codons}, amino_acids={amino_acids})'


DEFAULT_GENETIC_CODE = GeneticCode()

SUPE = GeneticCode(updated_mappings={'TAG': 'Q'}, name='supE')


def translate(dna, genetic_code=None):
    """
    Translates a DNA string into amino acid string.

    RNA and lower-case strings accepted.
    """
    assert isinstance(dna, str)

    # process string by replacing 'U's with 'T's and uppercasing
    dna = dna.upper().replace('U', 'T')

    if not isinstance(genetic_code, GeneticCode):
        genetic_code = DEFAULT_GENETIC_CODE

    codon_strings = [dna[3 * i: 3 * i + 3] for i in range(len(dna) // 3)]

    return ''.join([genetic_code[_] for _ in codon_strings])


def test():
    for k, v in SUPE.code.items():
        print(f"'{k}': '{v}',")


if __name__ == '__main__':
    test()
