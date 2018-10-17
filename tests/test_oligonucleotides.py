from dogma.nucleotides import (
    Nucleotide)

from dogma.oligonucleotides import (
    Oligonucleotide)

def test_oligonucleotide_initialization():
    A = Nucleotide('A')
    B = Nucleotide('B')
    o1 = Oligonucleotide([A,B]*12)
    assert o1.label == 'AB'*12

    o2 = Oligonucleotide('AAA')
