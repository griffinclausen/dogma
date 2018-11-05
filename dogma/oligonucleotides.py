import decimal
import math

import pandas as pd

from dogma import (
    DEFAULT_DECIMAL_PRECISION,
    GeneticCode,
    DEFAULT_GENETIC_CODE,
    translate,
    Nucleotide,
    Codon,
    combine_nucleotides
)


decimal.getcontext().prec = DEFAULT_DECIMAL_PRECISION


class Oligonucleotide:
    """
    A sequence of Nucleotides objects.
    5' --> 3' list
    """

    def __init__(self, data, genetic_code=None, auto_run=True):
        """
        Multiple input formats are acceptable
            Oligonucleotide('AAANNK')
            Oligonucleotide(list('AAANNK'))
            Oligonucleotide(list(map(Nucleotide, list('AAANNK'))))
            Oligonucleotide([{'A':1}, {'A':0.5, 'C':0.5}]
            Oligonucleotide(list(map(Nucleotide, list('AAANNK'))))
        """

        if not isinstance(genetic_code, GeneticCode):
            genetic_code = DEFAULT_GENETIC_CODE
        self.genetic_code = genetic_code

        if isinstance(data, (str, list)):
            if isinstance(data[0], Codon):
                self.bases = []
                for c in data:
                    self.bases += c.bases
            else:
                self.bases = [_
                              if isinstance(_, Nucleotide)
                              else Nucleotide(_)
                              for _ in data]
        elif isinstance(data, Codon):
            self.bases = data.bases
            self.genetic_code = data.genetic_code
        else:
            self.bases = None

        self.length = self.get_length()
        self.label = self.get_label()
        self.compact_label = self.get_compact_label()

        self.codon_strings = self.get_codon_strings()
        self.codons = self.get_codons()

        self.get_base_profile()
        self.get_codon_profile()
        self.get_amino_acid_profile()

        if auto_run:
            self.assess_degeneracy()

    def get_bases(self):
        """
        Returns list of Nucleotides objects
        """
        return self.bases

    def get_codons(self):
        """
        Returns list of Codon objects derived from oligonucleotide bases
        """
        base_triplets = [self.bases[3 * i: 3 * i + 3]
                         for i in range(self.length // 3)]
        return [Codon(_, self.genetic_code) for _ in base_triplets]

    def get_length(self):
        """
        Returns integer length of oligonucleotide, number of bases
        """
        return len(self.bases)

    def get_label(self):
        """
        Returns a string of nucleotide labels.
        """
        return ''.join(_.label for _ in self.bases)

    def get_codon_strings(self):
        """
        Returns list of codon strings.
        """
        return [self.label[3 * i: 3 * i + 3] for i in range(self.length // 3)]

    def get_compact_label(self):
        """
        Returns shorthand label
        Attempts to reduce identical codon labels with numbers
        Ex: NNKNNKNNKAAANNK --> NNK<3>AAANNK
        """
        count = 1
        previous_string = ''
        output = []
        for c in self.get_codon_strings():
            if c == previous_string:
                count += 1
            else:
                if count > 1:
                    output.append(f'<{count}>')
                previous_string = c
                output.append(c)
                count = 1
        if count > 1:
            output.append(f'<{count}>')
        return ''.join(output)

    def translate(self, genetic_code=None):
        """
        Translates oligonucleotide into protein.

        Parameters:
        genetic_code: GeneticCode object
        """
        if not isinstance(genetic_code, GeneticCode):
            genetic_code = self.genetic_code
        return ''.join([genetic_code[_] for _ in self.get_codon_strings()])

    def samples(self, k=1, output='oligonucleotide', by_codon=True):
        return [self.sample(output=output, by_codon=by_codon)
                for _ in range(k)]

    def sample(self, output='oligonucleotide', by_codon=True):
        """
        Returns a single non-degenerate oligonucleotide string.

        If by_codon is True, sampling is based on the codons attribute,
        otherwise sampling is performed on each base individually.

        """
        if by_codon:
            oligonucleotide_string = ''.join(c.sample() for c in self.codons)
        else:
            oligonucleotide_string = ''.join([b.sample() for b in self.bases])

        if output == 'protein':
            return translate(oligonucleotide_string, self.genetic_code)
        elif output == 'both':
            return oligonucleotide_string, translate(oligonucleotide_string,
                                                     self.genetic_code)
        else:
            return oligonucleotide_string

    def assess_degeneracy(self):
        """
        Performs the sequential and potentially computationally and memory
        intensive calculations related to generating the protein abundance
        profile.
        """

        self.get_amino_acid_degeneracy_profile()
        self.get_protein_degeneracy_table()
        self.get_degeneracy_table()
        self.get_size()

    def get_base_profile(self):
        """
        Generates and returns base profile

        Usage
        -----
        >> nnk = Codon('NNK')
        >> oligo = Oligonucleotide(nnk)
        >> oligo.get_base_profile()
            {0: {'A':1, 'C':1 , 'G':1. 'T':1},
             1: {'A':1, 'C':1 , 'G':1. 'T':1},
             2: {               'G':1. 'T':1}}
        """
        self.base_profile = [b.composition for b in self.bases]
        return self.base_profile

    def get_codon_profile(self):
        """
        Generates and returns codon profile

        Usage
        -----
        >> nnk = Codon('NNK')
        >> oligo = Oligonucleotide(nnk)
        >> oligo.get_codon_profile()
            {0: {'AAA':1, 'AAC':1, AAG':1, 'AAT':1, 'ACA':1, ...},
             1: {'AAA':1, 'AAC':1, AAG':1, 'AAT':1, 'ACA':1, ...},
             2: {                 'AAG':1, 'AAT':1,          ...}}

        """
        self.codon_profile = [c.composition for c in self.codons]
        return self.codon_profile

    def get_amino_acid_profile(self):
        """
        Generates and returns amino acid profile.

        List of dictionaries mapping amino acid to degeneracy (redundance)
        Indicies of list reflects position within protein

        Usage
        -----
        >> supE = GeneticCode(1, {'TAG': 'Q'})
        >> oligo = Oligonucleotide('NNKGCC', supE)
        >> oligo.get_amino_acid_profile()
            [{'A':2, 'C':1, D':1, 'E':1, 'F':1, ...},
             {'A':1                                }]
        """
        self.amino_acid_profile = [c.translate(dtype='decimal')
                                   for c in self.codons]
        return self.amino_acid_profile

    def get_amino_acid_degeneracy_profile(self):
        """
        Generates and returns amino acid degeneracy profile

        Usage
        -----
        >> supE = GeneticCode(1, {'TAG': 'Q'})
        >> oligo = Oligonucleotide('NNKAAA', supE)
        >> oligo.get_amino_acid_degeneracy_profile()
            [{1:11, 2:6 , 3:3},
             {1:1            }]

        Explanation of profile (for 'NNKAAA' oligo with supE translation)
        ----------------------
        The 0 index refers to the first amino acid (or codon) position in oligo
            This NNK codon encodes 20 (11 + 6 + 3) amino acids
                11 of these amino acids (CDEFHIKMNQ) are encoded by 1 codon
                 6                      (AGPQTV)                    2 codons
                 3                      (LRS)                       3 codons
        The 1 index refers to the next position, where 'AAA' encodes for 'K'
                 1 amino acid           (K)                         1 codon
        """

        self.amino_acid_degeneracy_profile = []
        for d in self.amino_acid_profile:
            self.amino_acid_degeneracy_profile.append(
                {v: list(d.values()).count(v) for v in set(d.values())})

        return self.amino_acid_degeneracy_profile

    def get_protein_degeneracy_table(self):
        """
        Generates and returns protein degeneracy profile

        Protein degeneracy is a dictionary where keys are all possible
        degeneracies of proteins derived from the degenerate oligonucleotide,
        and where values are the number of proteins in this group.

        Dictionary values reflect number of proteins within Oligonucleotid
        with a given degeneracy.

        Usage
        -----
        >> supE = GeneticCode(1, {'TAG': 'Q'})
        >> oligo = Oligonucleotide('NNK'*2, supE)
        >> oligo.get_protein_degeneracy_table()
            {1: 121,  # key=1*1, value=11*11
             2: 132,  # key=1*2, value=11*6*2 (*2 because key=1*2 and 2*1)
             3: 66,
             4: 36,
             6: 36,
             9: 9}

        Explanation of profile (for 'NNKNNK' oligo with supE translation)
        ----------------------
        121 proteins (CC, CD, CE, ...) have a degeneracy of 1
        131 proteins (AC, AD, AE, ...) have a degeneracy of 2
        ...
        There are 20*20 = 400 = (121 + 132 + 66 + 36 + 36 + 9) unique protein
        members of this oligonucleotide ensemble.
        Multiplying the degeneracy by the protein counts give the number of
        unique DNA counts
        """

        # for each df, create new dfs for each row of next df,
        # then concat, repeat
        dfs = [pd.DataFrame({'Degeneracy': list(_.keys()),
                             'Proteins': list(map(decimal.Decimal,
                                                  _.values()))},
                            dtype='object')
               for _ in self.amino_acid_degeneracy_profile]

        if not dfs:
            return pd.DataFrame()

        df = dfs[0]
        for next_df in dfs[1:]:
            df = pd.concat([df * next_df.iloc[i] for i in range(len(next_df))])
            df = df.groupby(['Degeneracy'], sort=False).sum().reset_index()

        self.protein_degeneracy_table = df.sort_values(by=['Degeneracy']) \
                                          .reset_index(drop=True)
        return self.protein_degeneracy_table

    def get_degeneracy_table(self):
        """
        Expands protein degeneracy table
        """
        df = self.protein_degeneracy_table.copy()
        df['Oligonucleotides'] = df.Degeneracy * df.Proteins

        df['DNA_Quantile'] = \
            df.Oligonucleotides.cumsum() / sum(df.Oligonucleotides)
        df['Protein_Quantile'] = \
            df.Proteins.cumsum() / sum(df.Proteins)
        self.degeneracy_table = df

        # defining additional aliases
        self.protein_quantiles = df['Protein_Quantile'].tolist()
        self.dna_quantiles = df['DNA_Quantile'].tolist()
        self.degeneracies = df['Degeneracy'].tolist()

        self.df = df

    def get_size(self, output='dna'):
        """
        Calculates and returns size of oligonucleotide ensemble
        """

        self.size_oligonucleotides = self.df['Oligonucleotides'].sum()
        self.size_proteins = self.df['Proteins'].sum()
        self.size = self.size_oligonucleotides

        if output.startswith('d'):
            return self.size_oligonucleotides
        elif output.startswith('p'):
            return self.size_proteins
        else:
            return (self.size_oligonucleotides, self.size_proteins)

    def get_gini_index(self):
        """
        Calculates Gini Index of protein ensemble.
        """
        x = [0] + self.protein_quantiles
        y = [0] + self.dna_quantiles
        n = len(self.protein_quantiles)

        gini = decimal.Decimal(1 / 2)
        for i in range(n):
            gini -= (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2

        self.gini_index = 2 * gini
        return self.gini_index

    def get_makowski_diversity(self):
        """
        **INCORRECT**
        Calculates diversity as define by Makowski and Soares
        https://doi.org/10.1093/bioinformatics/btg013
        d = 1/(N*SUM(Pi^2))
        """
        d = decimal.Decimal(0)
        P = self.size_proteins

        degeneracy = self.degeneracy_table['Degeneracy']
        protein_counts = self.degeneracy_table['Proteins']

        for x, y in zip(protein_counts, degeneracy):
            p = y / self.size_oligonucleotides
            d += x * (p ** 2)

        self.makowski_diversity = 1 / (P * d)
        return self.makowski_diversity

    def _get_alternative_makowski_diversity(self):
        """
        EQUIVALENT ANSWER AS self.get_makowski_diversity()
        Calculates diversity as define by Makowski and Soares
        https://doi.org/10.1093/bioinformatics/btg013
        Equation (2)
        d = 1/(N*PRODUCTi{SUMj{Pij^2}})
        Where i is position protein,
              j is each amino acids,
              Pij is prob of amino acid j at i
        """
        d = decimal.Decimal(1)

        for aa_dict in self.amino_acid_profile:
            tot = sum(aa_dict.values())
            d *= sum([(aa / tot)**2 for aa in aa_dict.values()])

        self._makowski_diversity = 1 / (self.size_proteins * d)
        return self._makowski_diversity

    def get_entropy(self, base=None):
        """
        Calculates entropy of abundance profile
        s = -SUM(pi*ln(pi))
        """

        degeneracy = self.degeneracy_table['Degeneracy']
        protein_counts = self.degeneracy_table['Proteins']

        print(degeneracy)
        print(protein_counts)
        return 0
        entropy = decimal.Decimal(0)

        # for each degeneracy group
        for x, y in zip(protein_counts, degeneracy):
            p = y / self.size_proteins
            entropy += x * p * decimal.Decimal(math.log(p))

            if base is not None:
                entropy /= math.log(base)

        self.entropy = -entropy
        return self.entropy

    def __repr__(self):
        return f'Oligonucleotide({self.label})'


def combine_oligonucleotides(*oligos, proportions=None, offsets=None):
    """
    Merges oligonucleotides by combining nucleotide compositions.
    """
    n = len(oligos)
    m = max([_.length for _ in oligos])

    if proportions is None:
        proportions = [1] * n
    else:
        assert len(proportions) == n
        assert all([isinstance(_, (int, float)) for _ in proportions])
        assert all([0 <= _ for _ in proportions])

    if offsets is None:
        offsets = [0] * n
    else:
        assert len(offsets) == len(oligos)
        assert all([isinstance(_, int) for _ in offsets])
        assert all([0 <= _ <= m for _ in offsets])

    nucs = []
    for i in range(m):  # for each nucleotide position
        ns = []
        ps = []
        for oligo, offset, proportion in zip(oligos, offsets, proportions):
            if i >= offset and i < (offset + oligo.length):
                ns.append(oligo.bases[i - offset])
                ps.append(proportion)
        nucs.append(combine_nucleotides(*ns, proportions=ps))
    return Oligonucleotide(nucs)


def reverse_complement(oligo):
    """
    """
    if isinstance(oligo, Oligonucleotide):
        oligo = oligo.label

    # TODO!


def get_nnk(n=1):
    """
    Simple helper function to easily generate NNK-based oligos.
    """
    supE = GeneticCode(1, {'TAG': 'Q'})
    oligo = Oligonucleotide('NNK' * n, supE)
    return oligo
