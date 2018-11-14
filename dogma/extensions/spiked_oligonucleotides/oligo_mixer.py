"""
Case 1:
    Want to randomize a residue in a protein.
        if want specific AA composition: NNN, NNK, TriNuc, ... NxNyNz, swiftlib
        if want certain mutation rate: USE THIS MODULE


Usage
-----
from dogma import (
    GeneticCode,
    Oligonucleotide)
from spiked import Spiker
"""
import os
from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt

from dogma import (
    Oligonucleotide,
    GeneticCode,
    Nucleotide,
    rescale,
    STANDARD_NUCLEOTIDES)

from dogma.extensions.spiked_oligonucleotides.utils import (
    prod,
    hamming_distance,
    manhattan_distance)


here = os.path.abspath(os.path.dirname(__file__))


class MixedOligo(Oligonucleotide):
    """
    An Oligonucleotide generated from the mixture of a template oligo and a
    spiked oligo pattern, combined at a particular proportion.
    """

    def __init__(self,
                 template,
                 spike_pattern,
                 proportions=0.5):
        self.template = template
        self.spike_pattern = spike_pattern

        if not isinstance(proportions, list):
            proportions = [proportions] * self.template.length
        self.proportions = proportions

        o1 = self.template
        o2 = self.spike_pattern
        ps = self.proportions

        bases = [combine_bases(i, j, proportions=[1-p, p])
                 for i, j, p in zip(o1.bases, o2.bases, ps)]

        super().__init__(bases, template.genetic_code)

        # self.calculate_mutation_rates()
        # self.bootstrap_dna_mutation_rate_count_profile()

    def calculate_mutation_rates(self):
        self.calculate_dna_mutation_rates()
        self.calculate_protein_mutation_rates()
        # self.calculate_truncation_rates()

    def calculate_dna_mutation_rates(self):
        """
        List of Manhattan (L1) distances between each aligned base composition,
        for each mixed oligonucleotide in self.mixed_oligos
        """
        xs = [_.composition for _ in self.template.bases]
        ys = [_.composition for _ in self.bases]

        self.dna_mutation_rates = [manhattan_distance(x, y)
                                   for x, y in zip(xs, ys)]
        self.calculate_dna_mutation_rate()
        return self.dna_mutation_rates

    def calculate_dna_mutation_rate(self):
        """
        Calculates overall rate of mutation, equal to the probability that
        a randomly selected oligo from the degenerate oligo design has at least
        one DNA mutation relative to the template oligo.
        """
        self.dna_mutation_rate = prod([1 - _ for _ in self.dna_mutation_rates])
        return self.dna_mutation_rate

    def calculate_dna_mutation_rate_count_profile(self):
        """
        Calculates the probability of a randomly sampled oligo from a the
        degenerate mixed_oligo will contain n mutations (where n in [0, L] and
        L is the length of the DNA oligo)

        TODO: may be too computationally complex
        """
        pass

    def bootstrap_dna_mutation_rate_count_profile(self, N=1000):
        """
        Samples N oligos from MixedOligo and groups by DNA mutation count
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

    def calculate_protein_mutation_rates(self):
        """
        List of Manhattan (L1) distances between each aligned residue
        composition for each mixed oligonucleotide in self.mixed_oligos
        """
        aa0 = self.amino_acid_profile
        aa1 = self.template.amino_acid_profile
        self.protein_mutation_rates = [manhattan_distance(x, y)
                                       for x, y in zip(aa0, aa1)]
        self.calculate_protein_mutation_rate()
        return self.protein_mutation_rates

    def calculate_protein_mutation_rate(self):
        """
        Calculates overall rate of mutation, equal to the probability that
        a randomly selected oligo from the degenerate oligo design has at least
        one protein mutation relative to the template oligo.
        """
        self.protein_mutation_rate = prod([1 - _ for _ in self.protein_mutation_rates])
        return self.protein_mutation_rate

    def bootstrap_protein_mutation_rate_count_profile(self, N=1000):
        """
        Samples N oligos from MixedOligo and groups by Protein mutation count
        """
        proteins = self.samples(N, output='protein')
        mutation_counts = [hamming_distance(_, self.template.translate())
                           for _ in proteins]
        mutation_profile = Counter(mutation_counts)
        return mutation_profile


class OligoMixer:
    """
    Generates new degenerates oligonucleotide from template oligo and an added
    oligonucleotide.
    """
    def __init__(self, template, added_oligo):
        self.template = template
        self.added_oligo = added_oligo

    def sweep(self, start=0, stop=101, step=1, scale=100):
        """
        Creates a list of mixed Oligonucleotides generated from the combination
        of a range of proportions (default ratios 0 thtough 100% by 1%)
        """
        self.ratios = [_ / scale for _ in range(start, stop, step)]
        self.make_mixed_oligos(self.ratios)

    def make_mixed_oligos(self, ratios):
        assert min(ratios) >= 0
        assert max(ratios) <= 1

        self.ratios = ratios
        self.mixed_oligos = [MixedOligo(self.template, self.added_oligo, ratio)
                             for ratio in ratios]

        return self.mixed_oligos


def combine_bases(b1, b2, proportions):
    """
    Adds two bases by weighting nucleotide compositions
    """
    data = {}
    for n in STANDARD_NUCLEOTIDES:
        data[n] = 0
        for nuc, p in zip([b1, b2], proportions):
            data[n] += rescale(nuc.composition).get(n, 0) * p
    return Nucleotide(rescale(data))


def XXX_to_YYY_sweep_protein_mutation_profile_plot():
    gc = GeneticCode(1, {'TAG': 'Q'})
    XXX = Oligonucleotide('AAA'*4, gc)
    YYY = Oligonucleotide('NNK'*4)
    mixer = OligoMixer(XXX, YYY)
    mixer.sweep(start=0, stop=101, step=2, scale=100)

    data = []
    for i, mo in enumerate(mixer.mixed_oligos):
        print(f'{i} of {len(mixer.mixed_oligos)}')
        mo.calculate_mutation_rates()
        mp = mo.bootstrap_protein_mutation_rate_count_profile(10000)
        _data = {_: mp[_] / sum(mp.values()) for _ in mp}
        _data['p'] = mo.proportions[0]
        data.append(_data)

    df = pd.DataFrame(data).set_index('p')
    # df.to_clipboard()
    df.plot(linestyle='', marker='o', markersize=2)
    plt.savefig(os.path.join(here, 'mutation_profile.png'))
    with open(os.path.join(here, 'caption.txt'), 'w') as f:
        f.write('Bootstrapped Mutation Profile\n')
        f.write('Expected amino acid composition of oligonucleotides for ')
        f.write(f'increasing proportions of {XXX.label} ')
        f.write(f'spiked into {YYY.label}')


def XXX_to_YYY_sweep_mutation_profile_plot():
    # gc = GeneticCode()
    gc = GeneticCode(1, {'TAG': 'Q'})
    XXX = Oligonucleotide('AAA', gc)
    YYY = Oligonucleotide('NNK')
    mixer = OligoMixer(XXX, YYY)
    mixer.sweep()

    data = []
    for i, mo in enumerate(mixer.mixed_oligos):
        mo.calculate_mutation_rates()
        mp = mo.bootstrap_dna_mutation_rate_count_profile(1000)
        _data = {_: mp[_] / sum(mp.values()) for _ in mp}
        _data['p'] = mo.proportions[0]
        data.append(_data)

    df = pd.DataFrame(data).set_index('p')
    df.plot(linestyle='', marker='o', markersize=2)
    plt.savefig(os.path.join(here, 'mutation_rates.png'))
    with open(os.path.join(here, 'caption.txt'), 'w') as f:
        f.write('Bootstrapped Mutation Profile\n')
        f.write('Expected number of mutations per oligonucleotide for ')
        f.write(f'increasing proportions of {XXX.label} ')
        f.write(f'spiked into {YYY.label}')


def AAA_to_CCC_sweep_mutation_profile():
    AAA = Oligonucleotide('AAA')
    CCC = Oligonucleotide('CCC')
    for x in range(0, 101):
        x = x/100
        mo = MixedOligo(AAA, CCC, x)
        mo.calculate_mutation_rates()
        mp = mo.bootstrap_dna_mutation_rate_count_profile(10000)


def AAA_plus_CCC_MixOligo_sweep(xs):
    AAA = Oligonucleotide('AAA')
    CCC = Oligonucleotide('CCC')
    for x in xs:
        print(f'\n{x}:')
        mo = MixedOligo(AAA, CCC, x)
        mo.calculate_mutation_rates()
        mo.bootstrap_dna_mutation_rate_count_profile()


def AAA_plus_CCC_MixOligo(x=0.5):
    AAA = Oligonucleotide('AAA')
    CCC = Oligonucleotide('CCC')
    mo = MixedOligo(AAA, CCC, x)
    mo.calculate_mutation_rates()
    mo.bootstrap_dna_mutation_rate_count_profile()


def AAA_plus_NNN_Xpct(x=0.05):
    AAA = Oligonucleotide('AAA')
    NNN = Oligonucleotide('NNN')
    mixer = OligoMixer(AAA, NNN)
    o = mixer.mix(x)
    A = o.get_base_profile()[0]['A']
    aa = o.get_amino_acid_profile()[0]
    mutation_rate = 1 - aa.get('K', 0)
    print(f'{x}: A:{A:.2f} Mutation Rate: {mutation_rate:.2f}')


def AAA_plus_NNN_sweep():
    supE = GeneticCode(1, {'TAG': 'Q'})
    AAA = Oligonucleotide('AAA', supE)
    NNN = Oligonucleotide('CCC')
    mixer = OligoMixer(AAA, NNN)
    mixer.sweep()
    mixer.calculate_mutation_rates()


def analysis():
    # AAA_plus_NNN_sweep()
    # AAA_plusx_NNN_Xpct()
    # AAA_plus_CCC_MixOligo()
    # AAA_plus_CCC_MixOligo_sweep([0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
    # AAA_to_CCC_sweep_mutation_profile()
    # XXX_to_YYY_sweep_mutation_profile_plot()
    XXX_to_YYY_sweep_protein_mutation_profile_plot()


def main():
    analysis()


if __name__ == '__main__':
    main()
