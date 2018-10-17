"""
This module defines package-wide variables and helper functions.

References
----------
    Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences
    http://www.sbcs.qmul.ac.uk/iubmb/misc/naseq.html
"""

from itertools import product
from random import choice


DEFAULT_DECIMAL_PRECISION = 200  # Can be increased (ex: 1000) without significant speed issues, but unnecessary memory usage??

NBCI_GENETIC_CODE_NAMES = {
    1: 'Standard',
    2: 'Vertebrate Mitochondrial',
    3: 'Yeast Mitochondrial',
    4: 'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate; Mitochondrial; Mycoplasma; Spiroplasma',
    5: 'Invertebrate Mitochondrial',
    6: 'Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear',
    9: 'Echinoderm Mitochondrial; Flatworm Mitochondrial',
    10: 'Euplotid Nuclear',
    11: 'Bacterial, Archaeal and Plant Plastid',
    12: 'Alternative Yeast Nuclear',
    13: 'Ascidian Mitochondrial',
    14: 'Alternative Flatworm Mitochondrial',
    15: 'Blepharisma Macronuclear',
    16: 'Chlorophycean Mitochondrial',
    21: 'Trematode Mitochondrial',
    22: 'Scenedesmus obliquus Mitochondrial',
    23: 'Thraustochytrium Mitochondrial',
    24: 'Pterobranchia Mitochondrial',
    25: 'Candidate Division SR1 and Gracilibacteria',
    26: 'Pachysolen tannophilus Nuclear Code'
}

NCBI_NUCLEOTIDES = 'ACGT'

# list all 64 standard codons ['AAA', 'AAC', ... 'TTG', 'TTT']
NCBI_CODONS = [''.join(_) for _ in product(NCBI_NUCLEOTIDES, repeat=3)]

NCBI_AMINO_ACIDS = {
    1: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF',
    2: 'KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF',
    3: 'KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF',
    4: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF',
    5: 'KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF',
    6: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF',
    9: 'NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF',
    10: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF',
    11: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF',
    12: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF',
    13: 'KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF',
    14: 'NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF',
    15: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF',
    16: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF',
    21: 'NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF',
    22: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF',
    23: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF',
    24: 'KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF',
    25: 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF',
    26: 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
}

STANDARD_NUCLEOTIDES = 'ACGT'

DEGENERATE_NUCLEOTIDE_CODE = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': 'AG',
    'Y': 'CT',
    'S': 'CG',
    'W': 'AT',
    'K': 'GT',
    'M': 'AC',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG',
    'N': 'ACGT'
}

DEGENERATE_NUCLEOTIDES = sorted(DEGENERATE_NUCLEOTIDE_CODE.keys())  # includes standard nucleotides (ACGT)

DEGENERATE_NUCLEOTIDE_CODE_REVERSED = {v: k for k, v in DEGENERATE_NUCLEOTIDE_CODE.items()}

DEGENERATE_NUCLEOTIDE_CODE_COMPOSITION = {code: {n: int(n in _) for n in STANDARD_NUCLEOTIDES} for code, _ in DEGENERATE_NUCLEOTIDE_CODE.items()}

NUCLEOTIDE_BASE_PAIRS = dict(zip('ACGT', 'TGCA'))

DEGENERATE_NUCLEOTIDE_PAIRS = {k: DEGENERATE_NUCLEOTIDE_CODE_REVERSED[''.join(sorted([NUCLEOTIDE_BASE_PAIRS[_] for _ in v]))] for k, v in DEGENERATE_NUCLEOTIDE_CODE.items()}

DEFAULT_NUCLEOTIDE_LABEL = 'N'

DEGENERATE_NUCLEOTIDE_PAIRS_ = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'R': 'Y',
    'Y': 'R',
    'S': 'S',
    'W': 'W',
    'K': 'M',
    'M': 'K',
    'B': 'V',
    'D': 'H',
    'H': 'D',
    'V': 'B',
    'N': 'N'
}


DEFAULT_OLIGONUCLEOTIDE_LABEL = 'NNN'

STANDARD_CODONS = [''.join(_) for _ in product(STANDARD_NUCLEOTIDES, repeat=3)]

DEGENERATE_CODONS = [''.join(_) for _ in product(DEGENERATE_NUCLEOTIDES, repeat=3)]

DEFAULT_CODON_LABEL = 'NNN'


STANDARD_AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

DEFAULT_AMINO_ACID_LABEL = 'X'

STOP_LABEL = '*'


DEFAULT_RESIDUE_LABEL = 'X'


def rescale(data, total=1):
    """
    Rescales numerical values in lists and dictionaries to sum to specified total.
    """

    if isinstance(data, list):
        input_total = sum(data)
        assert input_total != 0, 'Error in dogma.utils.rescale(), input_total == 0'

        return [_ / input_total * total for _ in data]

    elif isinstance(data, dict):
        input_total = sum(data.values())
        assert input_total != 0, 'Error in dogma.utils.rescale(), input_total == 0'

        return {k: v / input_total * total for k, v in data.items()}


def get_frequency_dictionary(data):
    """
    Takes a string or list of strings and returns a dictionary of unique members and their abundance.
    """
    return {k: data.count(k) for k in set(data)}


def get_random_oligonucleotide(length=3, letters='ACGTRYSWKMBDHVN'):
    """
    Returns string of random degenerate nucleotides.
    Each of the 15 letters 'ACGTRYSWKMBDHVN' are equally likely.
    Alternatively, use letters='ACGT' for standard oligo strings
    """
    return ''.join(choice(list(letters)) for _ in range(length))


def test_get_random_oligonucleotide():
    print(get_random_oligonucleotide())
    print(get_random_oligonucleotide(12))
    print(get_random_oligonucleotide(12, 'ACGT'))


def tests():
    test_get_random_oligonucleotide()


def main():
    tests()


if __name__ == '__main__':
    main()
