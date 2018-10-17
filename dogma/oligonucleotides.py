import decimal 
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from dogma import (DEFAULT_DECIMAL_PRECISION,
                   GeneticCode,
                   DEFAULT_GENETIC_CODE,
                   translate,
                   Nucleotide,
                   Codon,
                   rescale,
                   STANDARD_NUCLEOTIDES,

                   get_random_oligonucleotide)


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
                self.bases = [_ if isinstance(_, Nucleotide) else Nucleotide(_) for _ in data]
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
        base_triplets = [self.bases[3*i:3*i+3] for i in range(self.length//3)]
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
        return [self.label[3*i:3*i+3] for i in range(self.length // 3)]


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
        return [self.sample(output=output, by_codon=by_codon) for _ in range(k)]


    def sample(self, output='oligonucleotide', by_codon=True):
        """
        Returns a single non-degenerate oligonucleotide string.

        If by_codon is True, sampling is based on the codons attribute, otherwise
        sampling is performed on each base individually.

        """
        if by_codon:
            oligonucleotide_string = ''.join(c.sample() for c in self.codons)
        else:
            oligonucleotide_string = ''.join([b.sample() for b in self.bases])

        if output == 'protein':
            return translate(oligonucleotide_string, self.genetic_code)
        elif output == 'both':
            return oligonucleotide_string, translate(oligonucleotide_string, self.genetic_code)
        else:
            return oligonucleotide_string


    def assess_degeneracy(self):
        """
        Performs the sequential and potentially computationally and memory intensive calculations
        related to generating the protein abundance profile.
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
        self.amino_acid_profile = [c.translate(dtype='decimal') for c in self.codons]
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
             self.amino_acid_degeneracy_profile.append( {v: list(d.values()).count(v) for v in set(d.values())} )

        return self.amino_acid_degeneracy_profile


    def get_protein_degeneracy_table(self):
        """
        Generates and returns protein degeneracy profile

        Protein degeneracy is a dictionary where keys are all possible degeneracies of proteins derived
        from the degenerate oligonucleotide, and where values are the number of proteins in this group.

        Dictionary values reflect number of proteins within Oligonucleotide with a given degeneracy.

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
        121 proteins (ex: CC, CD, CE, ...) have a degeneracy of 1, each encoded by 1 oligonucleotide
        131          (ex: AC, AD, AE, ...)                      2 ,                2
        ...
        There are 20*20 = 400 = (121 + 132 + 66 + 36 + 36 + 9) unique protein members of this oligonucleotide ensemble
        Multiplying the degeneracy by the protein counts give the number of unique DNA counts
        """

        # for each df, create new dfs for each row of next df, then concat, repeat


        dfs = [pd.DataFrame({'Degeneracy': list(_.keys()),
                             'Proteins': list(map(decimal.Decimal, _.values()))},
                             dtype='object') for _ in self.amino_acid_degeneracy_profile]

        if not dfs:
            return pd.DataFrame()

        df = dfs[0]
        for next_df in dfs[1:]:
            df = pd.concat([df * next_df.iloc[i] for i in range(len(next_df))])
            df = df.groupby(['Degeneracy'], sort=False).sum().reset_index()
        
        self.protein_degeneracy_table = df.sort_values(by=['Degeneracy']).reset_index(drop=True)
        return self.protein_degeneracy_table


    def get_degeneracy_table(self):
        """
        Expands protein degeneracy table
        """
        df = self.protein_degeneracy_table.copy()
        df['Oligonucleotides'] = df.Degeneracy * df.Proteins

        df['DNA_Quantile'] = df.Oligonucleotides.cumsum()/sum(df.Oligonucleotides)
        df['Protein_Quantile'] = df.Proteins.cumsum()/sum(df.Proteins)
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

        gini = decimal.Decimal(1/2)
        for i in range(n):
            gini -=  (x[i+1]-x[i]) * (y[i+1]+y[i])/2

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


    def heatmap(self, filepath=None):
        if filepath is None:
            filepath = 'heatmap.png'

        # prepare data
        a = self.amino_acid_profile
        a = [rescale({_: d.get(_, 0) for _ in self.genetic_code.amino_acids}) for i, d in enumerate(a)]
        
        b = self.base_profile
        b = [rescale({_: d.get(_, 0) for _ in STANDARD_NUCLEOTIDES}) for i, d in enumerate(b)]

        A = pd.DataFrame(a).T.astype(float)
        B = pd.DataFrame(b).T.astype(float)

        A = A.sort_index(ascending=False)
        B = B.sort_index(ascending=False)

        # plot heatmap
        fig = plt.figure(1)
        # gridspec.GridSpec(3,1)

        # Amino Acid
        ax0 = plt.subplot2grid((6, 1), (2, 0), rowspan=4)
        plt.sca(ax0)
        c = ax0.pcolor(A, edgecolors='k', linewidths=1, cmap='Greens', vmin=0)
        ax0.set_title('Amino Acid Profile')
        plt.yticks(np.arange(0.5, len(A.index), 1), A.index,
                   fontsize=8, horizontalalignment='center')
        plt.xticks(np.arange(0.5, len(A.columns), 1), A.columns,
                   fontsize=8, horizontalalignment='center')
        cbar = fig.colorbar(c, ax=ax0)
        cbar.ax.tick_params(labelsize=8) 

        # Base
        ax1 = plt.subplot2grid((6, 1), (0, 0), rowspan=2)
        plt.sca(ax1)
        c = ax1.pcolor(B, edgecolors='k', linewidths=1, cmap='Greens', vmin=0)
        ax1.set_title('Base Profile')
        plt.yticks(np.arange(0.5, len(B.index), 1), B.index,
                   fontsize=8, horizontalalignment='center')
        plt.xticks(np.arange(0.5, len(B.columns), 1), B.columns,
                   fontsize=8, horizontalalignment='center')
        cbar = fig.colorbar(c, ax=ax1)
        cbar.ax.tick_params(labelsize=8) 

        fig.tight_layout()
        fig.savefig(filepath, dpi=150)


    def bar(self, filepath=None, sort_aa=False):
        if filepath is None:
            filepath = 'bar.png'

        # prepare data
        N = self.length // 3
        a = self.amino_acid_profile
        a = [rescale({_: d.get(_, 0) for _ in self.genetic_code.amino_acids}) for i, d in enumerate(a)]
        
        b = self.base_profile
        b = [rescale({_: d.get(_, 0) for _ in STANDARD_NUCLEOTIDES}) for i, d in enumerate(b)]

        A = pd.DataFrame(a).T.astype(float).sort_index(ascending=False)
        B = pd.DataFrame(b).T.astype(float).sort_index(ascending=False)

        # create figure
        if N > 5:
            figsize = (6.4*N/5, 4.8)
        else:
            figsize = (6.4, 4.8)
        fig = plt.figure(figsize=figsize)
        
        # For each codon
        for n in range(N):

            # Base
            data = B.ix[:, [3*n, 3*n+1, 3*n+2]]

            ax0 = plt.subplot2grid((3, N), (0, n))
            c = ax0.pcolor(data, edgecolors='k', linewidths=1, cmap='Greens', vmin=0, vmax=1)

            codon_label = ''.join([_.label for _ in self.bases[3*n:3*n+3]])
            
            plt.yticks(np.arange(0.5, len(data.index), 1), data.index,
                       fontsize=8, horizontalalignment='center')
            plt.xticks(np.arange(0.5, len(data.columns), 1), data.columns,
                       fontsize=8, horizontalalignment='center')

            ax0.set_title(codon_label)
            if n == 0:
                ax0.set_ylabel('Base Profile')
                ax0.set_xlabel('Position', fontsize=8)

            # adding value labels to base profile
            for i in range(3):  # columns
                for j, d in enumerate(data.index):  # rows
                    v = data.iloc[j, i]
                    val = f'{100*v:{0}.{3}}'

                    if val.endswith('.0'):
                        val = val[:-2]
                    if data.iloc[j, i] == 1:
                        val = '100'

                    ax0.text(i + 0.5, j + 0.5, val, ha='center', va='center', fontsize=8)


            # Amino Acid
            data = A[n]
            if sort_aa:
                data = data.sort_values()

            ax1 = plt.subplot2grid((3, N), (1, n), rowspan=2)
            data.plot(kind='barh', ax=ax1, color='green', width=.85)

            plt.xticks(fontsize=8)
            plt.yticks(range(len(data.index)), data.index,
                       fontsize=8, horizontalalignment='center')
            if n == 0:
                ax1.set_ylabel('Amino Acid Profile')
                ax1.set_xlabel('Proportion', fontsize=8)

        fig.tight_layout()
        fig.savefig(filepath, dpi=150)


    def __repr__(self):
        return f'Oligonucleotide({self.label})'


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
    oligo = Oligonucleotide('NNK'*n, supE)
    return oligo    


def test_compact_label():
    o = Oligonucleotide('NNKNNKNNKAAANNK')
    print(o.get_compact_label())


def test_gini():
    supE = GeneticCode(1, {'TAG': 'Q'})
    for i in range(1, 13):
        o = Oligonucleotide('NNN'*i, supE)
        g = o.get_gini_index()
        print(g)


def test_max_nnk():
    for i in [1, 2, 3, 4, 5]:
        o = Oligonucleotide('NNK'*i*6)
        print(i)
    print(o.get_compact_label())
    print(len(o.protein_quantiles))


def test_nxnynz():
    supE = GeneticCode(1, {'TAG': 'Q'})
    data = [{'A':26, 'C':21, 'G':27, 'T':26},
            {'A':30, 'C':20, 'G':25, 'T':25},
            {'G':50, 'T':50}]

    import time
    for i in range(1, 3):
        tic = time.time()
        o = Oligonucleotide(data * i, supE)
        print(i, time.time()-tic)
        print(o.get_entropy())

    print(o.compact_label)
    print(len(o.degeneracies))


def test_makowski():
    supE = GeneticCode(1, {'TAG': 'Q'})
    for i in range(1, 13):
        o = Oligonucleotide('NNK' * i, supE)
        print(o.get_makowski_diversity())


def test_entropy():
    supE = GeneticCode(1, {'TAG': 'Q'})
    for i in range(1, 13):
        o = Oligonucleotide('NNN' * i, supE)
        print(o.get_entropy())


def test_heatmap():
    supE = GeneticCode(1, {'TAG': 'Q'})
    data = [
            {'A':26, 'C':21, 'G':27, 'T':26},
            {'A':30, 'C':20, 'G':25, 'T':25},
            {'G':50, 'T':50},

            {'A':40, 'C':10, 'G':0, 'T':26},
            {'A':4, 'C':10, 'G':1, 'T':25},
            {'G':50, 'T':50},

            {'A':26, 'C':0, 'G':0, 'T':26},
            {'A':30, 'C':0, 'G':0, 'T':25},
            {'G':50, 'A':23},
            ]
    o = Oligonucleotide(data, supE)
    o.heatmap()

def test_bar():
    supE = GeneticCode(1, {'TAG': 'Q'})
    data = [
            {'A':26, 'C':21, 'G':27, 'T':26},
            {'A':30, 'C':20, 'G':25, 'T':25},
            {'G':50, 'T':50},

            {'A':40, 'C':10, 'G':0, 'T':26},
            {'A':4, 'C':10, 'G':1, 'T':25},
            {'G':50, 'T':50},

            {'A':40, 'C':10, 'G':0, 'T':26},
            {'A':4, 'C':10, 'G':1, 'T':25},
            {'G':50, 'T':50},

            {'A':26, 'C':0, 'G':0, 'T':26},
            {'A':30, 'C':0, 'G':0, 'T':25},
            {'G':50, 'A':23},
            ]
    data = 'NNNNNBNNK'
    data = get_random_oligonucleotide(30)
    o = Oligonucleotide(data, supE)
    o.bar()  



if __name__ == '__main__':
    # test_get_random_oligonucleotide()
    # test_gini()
    # test_max_nnk()
    # test_nxnynz()
    # test_makowski()
    # test_entropy()
    # test_heatmap()
    test_bar()