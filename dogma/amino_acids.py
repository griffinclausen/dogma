#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dogma import (
    GeneticCode,
    DEFAULT_GENETIC_CODE,
    Codon,
)


class AminoAcid:
    """
    Amino acids are the building blocks of proteins.

    Attributes
    ++++++++++
    label - str - One letter string label denoting acid acid identity
    members - list - List of one letter symbols for nonzero amino acids
    proportions
    codon_labels
    codons
    composition

    Methods
    +++++++
    get_synonymous_codons
    is_degenerate
    is_equimolar

    """

    def __init__(self,
                 data=None,
                 genetic_code=None,
                 triplet_string_is_dna=True):

        if not isinstance(genetic_code, GeneticCode):
            genetic_code = DEFAULT_GENETIC_CODE
        self.genetic_code = genetic_code

        # input data is string of length 1, an amino acid letter
        # AminoAcid('A')
        if isinstance(data, str) and len(data) == 1:

            self.label = data
            self.members = [data]
            self.proportions = [1]
            self.codon_labels = self.get_synonymous_codons()
            self.codons = [self.get_synonymous_codons('Codon')]
            self.composition = dict(zip(self.members, self.proportions))

        # input data is a string of length 3, a codon label
        # AminoAcid('GGG')
        elif (isinstance(data, str) and
              len(data) == 3 and
              triplet_string_is_dna):

            self.codons = [Codon(data, genetic_code)]
            self.codon_labels = self.codons[0].members
            self.composition = self.codons[0].translate()
            self.members = sorted(self.composition.keys())
            self.proportions = [self.composition[_] for _ in self.members]
            self.label = genetic_code[data]

        # input data is a Codon object
        # AminoAcid(Codon('GGG'))
        elif isinstance(data, Codon):

            self.codons = [data]
            self.codon_labels = self.codons[0].members
            self.composition = self.codons[0].translate()
            self.members = sorted(self.composition.keys())
            self.proportions = [self.composition[_] for _ in self.members]
            self.label = genetic_code[data.label]

        # input data is dictionary of amino acids and their proportions
        elif isinstance(data, dict) and \
                all([_ in self.genetic_code.amino_acids for _ in data]) and \
                all([isinstance(_, (int, float)) for _ in data.values()]):
            self.codons = []
            self.codon_labels = ['']
            self.composition = data
            self.members = sorted(self.composition.keys())
            self.proportions = [self.composition[_] for _ in self.members]
            if not self.is_degenerate():
                self.label = f'({self.members})'
            else:
                self.label = f'{self.members[0]}'

        self.degenerate = self.is_degenerate()
        self.equimolar = self.is_equimolar()

    def get_synonymous_codons(self, output_type='string'):
        """
        Returns list of triplet strings or Codon instances.

        Each synonymous codon encodes a nonzero amino acid.
        """
        codon_strings = []
        for member in self.members:
            codon_strings += self.genetic_code.code_reversed[member]

        if output_type == 'string':
            return codon_strings
        elif output_type == 'Codon':
            return Codon(codon_strings, self.genetic_code)
        elif output_type == 'Codons':
            return [Codon(_, self.genetic_code) for _ in codon_strings]
        return None

    def is_degenerate(self):
        """
        Amino acid is degenerate if there are multiple amino acid members
        have nonzero composition values.
        """
        return len(self.members) > 0

    def is_equimolar(self):
        """
        Amino acid is equimolar if all nonzero amino acid members
        have equal abundance.
        """
        return max(self.proportions) == min(self.proportions)

    def __str__(self):

        print(self.label)
        print(self.composition)

    def __repr__(self):
        l, c = self.label, self.composition
        return f'AminoAcid(label="{l}", composition={c})'
