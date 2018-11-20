import os
from collections import Counter

from dogma import (
    Oligonucleotide,
    GeneticCode,
    Nucleotide,
    rescale,
    STANDARD_NUCLEOTIDES)

from dogma.extensions.spiked_oligonucleotides.utils import (
    prod,
    hamming_distance,
    manhattan_distance,
    combine_bases)


here = os.path.abspath(os.path.dirname(__file__))


class SpikedOligo(Oligonucleotide):
    """
    An Oligonucleotide generated from the mixture of a template oligo and a
    spiked oligo pattern, combined at a particular proportion.
    """

    def __init__(self,
                 template,
                 spike_pattern,
                 proportions=0.5,
                 auto_run=False):

        # Oligonucleotides
        self.template = template
        self.spike_pattern = spike_pattern

        # Proportion of template-to-spiked oligo at each base position
        if not isinstance(proportions, list):
            proportions = [proportions] * self.template.length
        self.proportions = proportions

        # Calculates nucleotide base composition for this resulting mixed oligo
        o1 = self.template
        o2 = self.spike_pattern
        ps = self.proportions
        bases = [combine_bases(i, j, proportions=[1-p, p])
                 for i, j, p in zip(o1.bases, o2.bases, ps)]
        super().__init__(bases, template.genetic_code, auto_run=auto_run)

        if auto_run:
            self.calculate_mutation_rates()

    def calculate_mutation_rates(self):
        self.calculate_dna_mutation_rates()
        self.calculate_protein_mutation_rates()
        self.calculate_truncation_rates()

    # CALCULATING DNA MUTATIONS PROBABILITIES
    def calculate_dna_mutation_rates(self):
        """
        List of Manhattan (L1) distances between each aligned base composition,
        for each mixed oligonucleotide in self.mixed_oligos
        """
        xs = [_.composition for _ in self.template.bases]
        ys = [_.composition for _ in self.bases]

        self.dna_mutation_rates = [manhattan_distance(x, y) / 2
                                   for x, y in zip(xs, ys)]
        self.calculate_overall_dna_mutation_rate()
        return self.dna_mutation_rates

    def calculate_overall_dna_mutation_rate(self):
        """
        Calculates overall rate of mutation, equal to the probability that
        a randomly selected oligo from the degenerate oligo design has at least
        one DNA mutation relative to the template oligo.
        """
        self.overall_dna_mutation_rate = prod([1 - _ for _ in self.dna_mutation_rates])
        return self.overall_dna_mutation_rate

    def calculate_dna_mutation_rate_count_profile(self):
        """
        Calculates the probability of a randomly sampled oligo from a the
        degenerate mixed_oligo will contain n mutations (where n in [0, L] and
        L is the length of the DNA oligo)
        """
        pass

    def bootstrap_dna_mutation_rate_count_profile(self, N=1000):
        """
        Samples N oligos from SpikedOligo and groups by DNA mutation count
        """
        oligos = self.samples(N)
        mutation_counts = [hamming_distance(_, self.template.label)
                           for _ in oligos]
        self.dna_mutation_profile = Counter(mutation_counts)

        self.calculate_dna_mutation_proportion_by_count()
        return self.dna_mutation_profile

    def calculate_dna_mutation_proportion_by_count(self):
        """
        """
        pass

    # CALCULATING Protein MUTATIONS PROBABILITIES
    def calculate_protein_mutation_rates(self):
        """
        List of Manhattan (L1) distances between each aligned residue
        composition for each mixed oligonucleotide in self.mixed_oligos
        """
        aa0 = self.amino_acid_profile
        aa1 = self.template.amino_acid_profile
        self.protein_mutation_rates = [manhattan_distance(x, y) / 2
                                       for x, y in zip(aa0, aa1)]
        self.calculate_overall_protein_mutation_rate()
        return self.protein_mutation_rates

    def calculate_overall_protein_mutation_rate(self):
        """
        Calculates overall rate of mutation, equal to the probability that
        a randomly selected oligo from the degenerate oligo design has at least
        one protein mutation relative to the template oligo.
        """
        self.overall_protein_mutation_rate = 1 - prod([1 - _ for _ in self.protein_mutation_rates])
        return self.overall_protein_mutation_rate

    def bootstrap_protein_mutation_rate_count_profile(self, N=1000):
        """
        Samples N oligos from SpikedOligo and groups by Protein mutation count
        """
        proteins = self.samples(N, output='protein')
        mutation_counts = [hamming_distance(_, self.template.translate())
                           for _ in proteins]
        mutation_profile = Counter(mutation_counts)
        return mutation_profile

    # CALCULATING TRUNCATION PROBABILITIES
    def calculate_truncation_rates(self):
        """
        Calculates the overall rate of truncated oligos, rate of stop codons
        at each codon position, and oligo length profile (from N'terminal)
        """
        # normalize amino acid profile to sum to 1 for each residue
        aas = [rescale(_) for _ in self.amino_acid_profile]
        self.stop_rates = [_.get('*', 0) for _ in aas]
        # Frequency of coding codon
        self.aa_rates = [1-_ for _ in self.stop_rates]

        # overall truncation rate is 1 minus the probability of non-stops at
        # each codon position in the mixed oligo
        self.overall_truncation_rate = 1 - prod([1 - _
                                                 for _ in self.stop_rates])

        # Frequency of translated proteins of lengths 0 to oligo.length//3
        self.protein_length_profile = {i: s * prod(self.aa_rates[:i])
                                       for i, s in enumerate(self.stop_rates)}
        self.protein_length_profile[self.length//3] = 1 - self.overall_truncation_rate
