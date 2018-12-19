from dogma.proteins import Protein
from dogma import (
    Oligonucleotide,
    Codon,
    AminoAcid)


def test_Protein_from_amino_acid_string():
    aa_string = 'ACDEFGHIKLMNPQRSTVWY'
    p = Protein(aa_string)
    assert p.label == aa_string


def test_Protein_from_nucleotide_string():
    dna_string = 'aaacccgggttt'
    p = Protein(dna_string, data_is_dna=True)
    assert p.label == 'LPGF'


def test_Protein_from_oligonucleotide_object():
    dna_string = 'aaacccgggttt'
    oligo = Oligonucleotide(dna_string)
    p = Protein(oligo)
    assert p.label == 'LPGF'


def test_Protein_from_list_of_codon_object():
    dna_string = 'aaa'
    c = Codon(dna_string)
    p = Protein(c)
    assert p.label == 'L'


def test_Protein_from_amino_acid_object():
    aa_string = 'A'
    aa = AminoAcid(aa_string)
    p = Protein(aa)
    assert p.label == aa_string


def test_Protein_from_amino_acid_list():
    aa_string = 'A'
    aa = [AminoAcid(aa_string)]
    p = Protein(aa)
    assert p.label == aa_string

    aa_string = 'ACDEFGHIKLMNPQRSTVWY'
    aas = [AminoAcid(_) for _ in aa_string]
    p = Protein(aas)
    assert p.label == aa_string


def test_Protein_from_codon_list():
    dna_string = 'aaa'
    c = [Codon(dna_string)]
    p = Protein(c)
    assert p.label == 'L'

    dna_string = 'aaacccgggttt'
    oligo = Oligonucleotide(dna_string)
    c = oligo.codons
    p = Protein(c)
    assert p.label == 'LPGF'
