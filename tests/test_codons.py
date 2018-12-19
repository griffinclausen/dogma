from dogma.codons import (
    Codon,
    combine_codons,
    nucleotide_string_to_codons,
    degenerate_codon_string_to_standard_members)
from dogma.nucleotides import Nucleotide
from dogma.oligonucleotides import Oligonucleotide


def test_Codon_string():

    c = Codon('AAA')
    assert c.label == 'AAA'
    assert len(c.composition) == 1
    assert c.members == ['AAA']
    assert c.proportions == [1]

    c = Codon('NNK')
    assert c.label == 'NNK'
    assert len(c.composition) == 32
    assert len(c.members) == 32
    assert c.proportions.count(1) == 32


def test_Codon_list_of_strings():

    c = Codon(list('AAA'))
    assert c.label == 'AAA'
    assert len(c.composition) == 1
    assert c.members == ['AAA']
    assert c.proportions == [1]

    c = Codon(list('NNK'))
    assert c.label == 'NNK'
    assert c.proportions.count(1) == 32


def test_Codon_list_of_codon_strings():

    c = Codon(['AAA'])
    assert len(c.composition) == 1
    assert c.members == ['AAA']
    assert c.proportions == [1]

    c = Codon(['AAA', 'CCC'])
    assert c.proportions.count(1) == 2


def test_Codon_list_of_nucleotides():

    o = Oligonucleotide('AAA')
    b = o.bases
    assert b[0].label == 'A'
    assert isinstance(b, list)
    assert isinstance(b[0], Nucleotide)
    assert isinstance(b[1], Nucleotide)
    assert isinstance(b[2], Nucleotide)

    c = Codon(b)
    assert c.label == 'AAA'
    assert len(c.composition) == 1
    assert c.proportions.count(1) == 1

    c = Codon(list(map(Nucleotide, 'NNK')))
    assert c.label == 'NNK'
    assert len(c.composition) == 32
    assert c.proportions.count(1) == 32


def test_combine_codons():
    c1 = Codon('AAA')
    c2 = Codon('AKN')
    c3 = combine_codons(c1, c2)
    assert c3.label == 'ADN'


def test_nucleotide_string_to_codons():
    s = 'AAAAAA'
    assert nucleotide_string_to_codons(s, 'strings') == ['AAA', 'AAA']


def test_degenerate_codon_string_to_standard_members():
    assert degenerate_codon_string_to_standard_members('NAA') == \
           ['AAA', 'CAA', 'GAA', 'TAA']
