import os

import pandas as pd
import matplotlib.pyplot as plt

from dogma import (
    Oligonucleotide,
    GeneticCode
)

from dogma.extensions.spiked_oligonucleotides.mixed_oligonucleotides import (
    SpikedOligo)

from dogma.extensions.spiked_oligonucleotides.oligo_mixer import (
    OligoMixer)


here = os.path.abspath(os.path.dirname(__file__))


def XXX_to_YYY_sweep_protein_mutation_bootstrapped_profile_plot():
    gc = GeneticCode(1, {'TAG': 'Q'})
    XXX = Oligonucleotide('AAA'*4, gc)
    YYY = Oligonucleotide('NNK'*4)
    mixer = OligoMixer(XXX, YYY)
    mixer.sweep(start=0, stop=101, step=2, scale=100)

    data = []
    for i, so in enumerate(mixer.spiked_oligos):
        print(f'{i} of {len(mixer.spiked_oligos)}')
        so.calculate_mutation_rates()
        mp = so.bootstrap_protein_mutation_rate_count_profile(10000)
        _data = {_: mp[_] / sum(mp.values()) for _ in mp}
        _data['p'] = so.proportions[0]
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


def XXX_to_YYY_sweep_protein_mutation_profile_plot():
    gc = GeneticCode(1, {'TAG': 'Q'})
    XXX = Oligonucleotide('AAACCCGGGTTT', gc)
    YYY = Oligonucleotide('NNK'*4)
    mixer = OligoMixer(XXX, YYY)
    mixer.sweep(start=0, stop=101, step=1, scale=100)

    data = []
    for i, so in enumerate(mixer.spiked_oligos):
        print(f'{i} of {len(mixer.spiked_oligos)}')
        # CHANGE BOOTSTRAP
        mp = so.bootstrap_protein_mutation_rate_count_profile(1000)
        _data = {_: mp[_] / sum(mp.values()) for _ in mp}
        _data['p'] = so.proportions[0]
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
    XXX = Oligonucleotide('AAA'*3, gc)
    YYY = Oligonucleotide('NNK'*3)
    mixer = OligoMixer(XXX, YYY)
    mixer.sweep()

    data = []
    for i, so in enumerate(mixer.spiked_oligos):
        so.calculate_mutation_rates()
        mp = so.bootstrap_dna_mutation_rate_count_profile(1000)
        _data = {_: mp[_] / sum(mp.values()) for _ in mp}
        _data['p'] = so.proportions[0]
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
        so = SpikedOligo(AAA, CCC, x)
        so.calculate_mutation_rates()
        mp = so.bootstrap_dna_mutation_rate_count_profile(10000)


def AAA_plus_CCC_MixOligo_sweep(xs):
    AAA = Oligonucleotide('AAA')
    CCC = Oligonucleotide('CCC')
    for x in xs:
        print(f'\n{x}:')
        so = SpikedOligo(AAA, CCC, x)
        so.calculate_mutation_rates()
        so.bootstrap_dna_mutation_rate_count_profile()


def AAA_plus_CCC_MixOligo(x=0.5):
    AAA = Oligonucleotide('AAA')
    CCC = Oligonucleotide('CCC')
    so = SpikedOligo(AAA, CCC, x)
    so.calculate_mutation_rates()
    so.bootstrap_dna_mutation_rate_count_profile()


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


def main():
    # AAA_plus_NNN_sweep()
    # AAA_plusx_NNN_Xpct()
    # AAA_plus_CCC_MixOligo()
    # AAA_plus_CCC_MixOligo_sweep([0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
    # AAA_to_CCC_sweep_mutation_profile()
    # XXX_to_YYY_sweep_mutation_profile_plot()
    # XXX_to_YYY_sweep_protein_mutation_bootstrapped_profile_plot()
    XXX_to_YYY_sweep_protein_mutation_profile_plot()


if __name__ == '__main__':
    main()
