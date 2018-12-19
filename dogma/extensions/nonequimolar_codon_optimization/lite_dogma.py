from itertools import product

from dogma import STANDARD_CODONS
from dogma.utils import (
    product_of_list,
    rescale
)


class LiteGeneticCode:
    """
    Light-weight GeneticCode for efficient optimization.
    """
    def __init__(self, genetic_code):
        self.amino_acids = sorted(list(set(genetic_code.amino_acids)))
        self.aa_codon_idxs = [[i for i, c in enumerate(STANDARD_CODONS)
                               if genetic_code[c] == aa]
                              for aa in self.amino_acids]


class LiteBase:
    """
    Light-weight Base class for efficient optimization.
    """
    def __init__(self, acgt):
        self.acgt = acgt  # [a, c, g, t]
        self.acgt = rescale(self.acgt, 100)


class LiteCodon:
    """
    Light-weight Codon class for efficient optimization.
    """
    def __init__(self, bases, genetic_code):
        self.bases = bases
        self.genetic_code = genetic_code

        self.amino_acids = genetic_code.amino_acids

        self.base_profile = self.get_base_profile()
        self.codon_profile = self.get_codon_profile()
        self.amino_acid_profile = self.get_amino_acid_profile()

    def get_base_profile(self):
        """
        list of lists: [[a1, c1, g1, t1], [a2, c2, g2, t2], [a3, c3, g3, t3]]
        """
        return [_.acgt for _ in self.bases]

    def get_codon_profile(self):
        """
        list of codon counts (in alphabetical order of all 64 STANDARD_CODONS)
        """
        return [product_of_list(_) for _ in product(*self.base_profile)]

    def get_amino_acid_profile(self):
        """
        list of amino acid counts
        (in alphabetical order for all amino acids in genetic code scope)
        """
        return [sum([self.codon_profile[i] for i in idxs])
                for idxs in self.genetic_code.aa_codon_idxs]
