from dogma.nucleotides import (
    Nucleotide)

from dogma.oligonucleotides import (
    Oligonucleotide,
    combine_oligonucleotides,
    reverse_complement,
    calculate_protein_degeneracy,
    calculate_protein_quantile)

from dogma import (
    GeneticCode
)


def test_oligonucleotide_initialization():
    A = Nucleotide('A')
    B = Nucleotide('B')
    o1 = Oligonucleotide([A, B]*12)
    assert o1.label == 'AB'*12

    o2 = Oligonucleotide('AAA')
    assert isinstance(o2, Oligonucleotide)


def test_combine_oligonucleotides():
    G = Oligonucleotide('G')
    T = Oligonucleotide('T')
    K = Oligonucleotide('K')

    k = combine_oligonucleotides(G, T)
    assert isinstance(k, Oligonucleotide)
    assert k.get_base_profile(rescaled=1) == K.get_base_profile(rescaled=1)


def test_reverse_complement():
    oligo_string = 'ACCGGGTTTT'
    o1 = Oligonucleotide(oligo_string)
    o2 = reverse_complement(o1)
    assert o2.label == 'AAAACCCGGT'


def test_calculate_protein_degeneracy():
    supE = GeneticCode(1, {'TAG': 'Q'})
    nnn = Oligonucleotide('NNN', supE)
    nnk = Oligonucleotide('NNK', supE)
    aaa = Oligonucleotide('AAA', supE)
    nnk7 = Oligonucleotide('NNK'*7, supE)

    assert calculate_protein_degeneracy('A', nnn) == 4
    assert calculate_protein_degeneracy('A', nnk) == 2
    assert calculate_protein_degeneracy('A', aaa) == 0

    assert calculate_protein_degeneracy('C'*7, nnk7) == 1
    assert calculate_protein_degeneracy('G'*7, nnk7) == 2**7
    assert calculate_protein_degeneracy('R'*7, nnk7) == 3**7


def test_calculate_protein_quantile():
    supE = GeneticCode(1, {'TAG': 'Q'})
    nnn = Oligonucleotide('NNN', supE)
    nnk = Oligonucleotide('NNK', supE)
    aaa = Oligonucleotide('AAA', supE)
    nnk7 = Oligonucleotide('NNK'*7, supE)

    assert 0.856 < calculate_protein_quantile('A', nnn) < 0.858
    assert 0.849 < calculate_protein_quantile('A', nnk) < 0.851
    assert calculate_protein_quantile('A', aaa) == -1

    assert 0.0151 < calculate_protein_quantile('C'*7, nnk7) < 0.0153
    assert 0.969 < calculate_protein_quantile('G'*7, nnk7) < 0.971
    assert .999 < calculate_protein_quantile('R'*7, nnk7) < 1.001
