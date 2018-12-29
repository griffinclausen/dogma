from dogma.amino_acids import AminoAcid

from dogma import Codon


def test_AminoAcid_from_letter():
    A = AminoAcid('A')
    assert A.label == 'A'
    assert A.members == ['A']
    assert A.proportions == [1]
    assert len(A.codons) == 1 and all([isinstance(_, Codon) for _ in A.codons])
    assert len(A.codons[0].members) == 4
    assert A.codon_labels == ['GCA', 'GCC', 'GCG', 'GCT']
    assert A.composition == {'A': 1}


def test_AminoAcid_from_codon_string():
    G = AminoAcid('GGG')
    assert G.label == 'G'
    assert G.members == ['G']
    assert G.proportions == [1]
    assert len(G.codons) == 1 and all([isinstance(_, Codon) for _ in G.codons])
    assert G.codon_labels == ['GGG']
    assert G.composition == {'G': 1}


def test_AminoAcid_from_standard_codon():
    ggg = Codon('GGG')
    G = AminoAcid(ggg)
    assert G.label == 'G'
    assert G.members == ['G']
    assert G.proportions == [1]
    assert len(G.codons) == 1 and all([isinstance(_, Codon) for _ in G.codons])
    assert G.codon_labels == ['GGG']
    assert G.composition == {'G': 1}


def test_AminoAcid_from_degenerate_codon():
    nnk = Codon('NNK')
    n = AminoAcid(nnk)
    assert n.label == '_'
    assert n.members == list('*ACDEFGHIKLMNPQRSTVWY')
    assert n.proportions == [1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3,
                             1, 1, 2, 1, 3, 3, 2, 2, 1, 1]
    assert len(n.codons) == 1 and all([isinstance(_, Codon) for _ in n.codons])
    assert n.composition == dict(zip('RLSGAPTVQCDEFHIKMNYW*',
                                     list(map(int, '333222221111111111111'))))
