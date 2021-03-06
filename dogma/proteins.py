#!/usr/bin/env python
# -*- coding: utf-8 -*-

from functools import reduce
from random import choices

from dogma import (
    GeneticCode,
    DEFAULT_GENETIC_CODE,
    Oligonucleotide,
    Codon,
    AminoAcid,
)


class Protein:
    """
    A sequence of AminoAcid objects.

    N-terminal --> C-terminal list

    Parameters
    ----------
    data:
    genetic_code: GeneticCode object
    data_is_dna: flag for processing data

    Attributes
    ---------
    oligonucleotide:
    codons
    residues
    label
    length

    Methods
    -------
    is_degenerate
    is_equimolar
    get_synonymous_codons
    sample
    samples
    """

    def __init__(self, data, genetic_code=None, data_is_dna=False):

        if not isinstance(genetic_code, GeneticCode):
            genetic_code = DEFAULT_GENETIC_CODE
        self.genetic_code = genetic_code


        # Protein(Oligonucleotide('NNK'))
        if isinstance(data, Oligonucleotide):
            self.oligonucleotide = data
            self.codons = self.oligonucleotide.codons
            self.residues = [AminoAcid(c.translate()) for c in self.codons]

        # Protein(Codon('NNK'))
        elif isinstance(data, Codon):
            self.oligonucleotide = Oligonucleotide(data.bases, data.genetic_code)
            self.codons = self.oligonucleotide.codons
            self.residues = [AminoAcid(c.translate()) for c in self.codons]

        # Protein(AminoAcid('A'))
        elif isinstance(data, AminoAcid):
            self.residues = [data]
            self.oligonucleotide = None
            self.codons = None

        elif all([isinstance(_, AminoAcid) for _ in data]):
            self.residues = data
            self.oligonucleotide = None
            self.codons = None

        elif all([isinstance(_, Codon) for _ in data]):
            self.codons = data
            self.residues = [AminoAcid(_, _.genetic_code) for _ in data]
            self.oligonucleotide = None

        elif all([isinstance(_, str) for _ in data]):

            # Protein('NNK', data_is_dna=True)
            if data_is_dna:
                self.oligonucleotide = Oligonucleotide(data, genetic_code)
                self.codons = self.oligonucleotide.codons
                self.residues = [AminoAcid(c.translate()) for c in self.codons]

            # Protein('NNK')
            else:
                self.residues = [AminoAcid(_, genetic_code) for _ in data]
                self.codons = [a.get_synonymous_codons(Codon)
                               for a in self.residues]
                self.oligonucleotide = None

        else:
            self.residues = []
            self.oligonucleotide = None
            self.codons = []

        self.label = ''.join([_.label for _ in self.residues])
        self.length = len(self.residues)
        self.degenerate = self.is_degenerate()

    def is_degenerate(self):
        return any([a.is_degenerate() for a in self.residues])

    def get_synonymous_codons(self):
        return [codon for codon in [a.get_synonymous_codons()
                for a in self.amino_acids]]

    def samples(self, k=1):
        return [self.sample() for _ in range(k)]

    def sample(self):
        return choices(self.letters, self.proportions)

    def __str__(self):
        return self.label

    def __repr__(self):
        return f'Protein(label={self.label})'


def calculate_protein_degeneracy(protein, oligonucleotide=None):
    """
    Calculates degeneracy of protein sequences for a specified
    oligonucleotide design.

    If no oligonucleotide is specified, defaults to a fully randomized
    ('NNN') scheme with a length 3 times as large as the protein, and
    the default genetic code is used.
    """

    if isinstance(protein, Protein):
        protein = protein.label  # how to handle '_'?

    if oligonucleotide is None:
        oligonucleotide = Oligonucleotide('NNN' * len(protein))

    data = oligonucleotide.amino_acid_profile

    return reduce(lambda x, y: x * y, [aa_composition[aa]
                  for aa_composition, aa in zip(data, protein)],)
