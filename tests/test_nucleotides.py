from dogma.nucleotides import(
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


def test_is_nondegenerate_nucleotide_string():
    
    assert is_nondegenerate_nucleotide_string('A') == True
    assert is_nondegenerate_nucleotide_string('a') == True
    assert is_nondegenerate_nucleotide_string('T') == True
    assert is_nondegenerate_nucleotide_string('U') == True
    assert is_nondegenerate_nucleotide_string('u') == True

    assert is_nondegenerate_nucleotide_string('N') == False
    assert is_nondegenerate_nucleotide_string('n') == False
    assert is_nondegenerate_nucleotide_string('K') == False
    
    assert is_nondegenerate_nucleotide_string('AAA') == True
    assert is_nondegenerate_nucleotide_string('ANN') == False



def test_is_degenerate_nucleotide_string():

    assert is_degenerate_nucleotide_string('A') == False
    assert is_degenerate_nucleotide_string('a') == False
    assert is_degenerate_nucleotide_string('T') == False
    assert is_degenerate_nucleotide_string('U') == False
    assert is_degenerate_nucleotide_string('u') == False

    assert is_degenerate_nucleotide_string('N') == True
    assert is_degenerate_nucleotide_string('n') == True
    assert is_degenerate_nucleotide_string('K') == True

    assert is_degenerate_nucleotide_string('AAA') == False
    assert is_degenerate_nucleotide_string('AAN') == True


def test_is_valid_nucleotide_string():

    s = 'ACGTACGTAGCTAGCTGACG'
    assert is_valid_nucleotide_string(s) == True

    s = 'NWNMBHV'
    assert is_valid_nucleotide_string(s) == True

    s = ''
    assert is_valid_nucleotide_string(s) == True

    s = 'AAAAAAZ'
    assert is_valid_nucleotide_string(s) == False

    s = 'acgtacgtagctagctgacg'
    assert is_valid_nucleotide_string(s) == True

    s = 'nwnmbhv'
    assert is_valid_nucleotide_string(s) == True


def test_combine():

    c = Nucleotide('C')
    t = Nucleotide('T')

    c_plus_t = combine_nucleotides(c, t)
    assert isinstance(c_plus_t, Nucleotide)
    assert c_plus_t.composition == {'A': 0, 'C': 0.5, 'G': 0, 'T':0.5}
    assert c_plus_t.members == 'CT'
    assert c_plus_t.label == 'Y'
    assert c_plus_t.is_degenerate() == True
    assert c_plus_t.is_equimolar() == True

    c_plus_ttt = combine_nucleotides(c, t, proportions=[1, 3])
    assert isinstance(c_plus_ttt, Nucleotide)
    assert c_plus_ttt.members == 'CT'
    assert c_plus_ttt.label == 'y'
    assert c_plus_ttt.is_degenerate() == True
    assert c_plus_ttt.is_equimolar() == False
    assert c_plus_ttt.composition == {'A': 0, 'C': 0.25, 'G': 0, 'T':0.75}


def test_nucleotide_composition_to_letter():
    d = {'A': 1}
    assert nucleotide_composition_to_letter(d) == 'A' 

    d = dict(zip('ACGT', [1, 1, 1, 1]))
    assert nucleotide_composition_to_letter(d) == 'N' 

    d = dict(zip('ACGT', [1, 1, 2, 1]))
    assert nucleotide_composition_to_letter(d) == 'n' 
