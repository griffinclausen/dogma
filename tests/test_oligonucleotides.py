from dogma.nucleotides import (
    Nucleotide)

from dogma.oligonucleotides import (
    Oligonucleotide,
    combine_oligonucleotides)


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
