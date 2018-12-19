from dogma.nucleotides import (
    Nucleotide,
    is_nondegenerate_nucleotide_string,
    is_degenerate_nucleotide_string,
    is_valid_nucleotide_string,
    nucleotide_composition_to_letter,
    combine_nucleotides,
)


def test_nucleotide_initialization():

    A = Nucleotide('A')
    C = Nucleotide('C')
    N = Nucleotide('N')

    assert A.label == 'A'
    assert C.label == 'C'
    assert N.label == 'N'


def test_is_nondegenerate_nucleotide_string():

    assert is_nondegenerate_nucleotide_string('A') is True
    assert is_nondegenerate_nucleotide_string('a') is True
    assert is_nondegenerate_nucleotide_string('T') is True
    assert is_nondegenerate_nucleotide_string('U') is True
    assert is_nondegenerate_nucleotide_string('u') is True

    assert is_nondegenerate_nucleotide_string('N') is False
    assert is_nondegenerate_nucleotide_string('n') is False
    assert is_nondegenerate_nucleotide_string('K') is False

    assert is_nondegenerate_nucleotide_string('AAA') is True
    assert is_nondegenerate_nucleotide_string('ANN') is False


def test_is_degenerate_nucleotide_string():

    assert is_degenerate_nucleotide_string('A') is False
    assert is_degenerate_nucleotide_string('a') is False
    assert is_degenerate_nucleotide_string('T') is False
    assert is_degenerate_nucleotide_string('U') is False
    assert is_degenerate_nucleotide_string('u') is False

    assert is_degenerate_nucleotide_string('N') is True
    assert is_degenerate_nucleotide_string('n') is True
    assert is_degenerate_nucleotide_string('K') is True

    assert is_degenerate_nucleotide_string('AAA') is False
    assert is_degenerate_nucleotide_string('AAN') is True


def test_is_valid_nucleotide_string():

    s = 'ACGTACGTAGCTAGCTGACG'
    assert is_valid_nucleotide_string(s) is True

    s = 'NWNMBHV'
    assert is_valid_nucleotide_string(s) is True

    s = ''
    assert is_valid_nucleotide_string(s) is True

    s = 'AAAAAAZ'
    assert is_valid_nucleotide_string(s) is False

    s = 'acgtacgtagctagctgacg'
    assert is_valid_nucleotide_string(s) is True

    s = 'nwnmbhv'
    assert is_valid_nucleotide_string(s) is True


def test_combine():

    c = Nucleotide('C')
    t = Nucleotide('T')

    c_plus_t = combine_nucleotides(c, t)
    assert isinstance(c_plus_t, Nucleotide)
    assert c_plus_t.composition == {'A': 0, 'C': 0.5, 'G': 0, 'T': 0.5}
    assert c_plus_t.members == 'CT'
    assert c_plus_t.label == 'Y'
    assert c_plus_t.is_degenerate() is True
    assert c_plus_t.is_equimolar() is True

    c_plus_ttt = combine_nucleotides(c, t, proportions=[1, 3])
    assert isinstance(c_plus_ttt, Nucleotide)
    assert c_plus_ttt.members == 'CT'
    assert c_plus_ttt.label == 'y'
    assert c_plus_ttt.is_degenerate() is True
    assert c_plus_ttt.is_equimolar() is False
    assert c_plus_ttt.composition == {'A': 0, 'C': 0.25, 'G': 0, 'T': 0.75}


def test_nucleotide_composition_to_letter():
    d = {'A': 1}
    assert nucleotide_composition_to_letter(d) == 'A'

    d = dict(zip('ACGT', [1, 1, 1, 1]))
    assert nucleotide_composition_to_letter(d) == 'N'

    d = dict(zip('ACGT', [1, 1, 2, 1]))
    assert nucleotide_composition_to_letter(d) == 'n'
