#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Translation system information from NCBI.

:Reference:
    Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences
    http://www.sbcs.qmul.ac.uk/iubmb/misc/naseq.html
"""

from itertools import product


# dictionary of NCBI genetic code ID numbers and full names
NBCI_GENETIC_CODE_NAMES = {
    1: 'Standard',
    2: 'Vertebrate Mitochondrial',
    3: 'Yeast Mitochondrial',
    4: 'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate;'
       'Mitochondrial; Mycoplasma; Spiroplasma',
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

# string of one-letter codes for the four standard nucleotides
NCBI_NUCLEOTIDES = 'ACGT'

# list of all 64 standard codons ['AAA', 'AAC', ... 'TTG', 'TTT']
NCBI_CODONS = [''.join(_) for _ in product(NCBI_NUCLEOTIDES, repeat=3)]

# string of one-letter amino acid codes aligned with alphabetized codon list
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
