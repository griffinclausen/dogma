import pandas as pd
from dogma import (
    DEGENERATE_CODONS
)

from dogma.utils import (
    DEGENERATE_NUCLEOTIDE_CODE,
    rescale
)

from dogma.extensions.nonequimolar_codon_optimization.lite_dogma import (
    LiteCodon,
    LiteBase
)

from dogma.extensions.nonequimolar_codon_optimization.config import (
    DEFAULTS
)
from dogma.extensions.nonequimolar_codon_optimization.utils import (
    get_random_base_profile
)


class NonEquimolarOptimizer:
    def __init__(self, **kwargs):

        # unpack user-defined parameters
        self.unpack_parameters(kwargs)

        # determine_optimal_equimolar_degenerate_codon
        self.best_dcs = self.find_best_equimolar_degenerate_codon()

    def unpack_parameters(self, params):
        self.genetic_code = params.get('genetic_code',
                                       DEFAULTS['genetic_code'])

        self.amino_acid_targets = params.get('amino_acid_targets',
                                             DEFAULTS['amino_acid_targets'])

        self.optimization_mode = params.get('optimization_mode',
                                            DEFAULTS['optimization_mode'])

        self.max_total_samples = params.get('max_total_samples',
                                            DEFAULTS['max_total_samples'])

        # create ordered list of amino acid targets based on genetic code scope
        self.targets = rescale([self.amino_acid_targets.get(_, 0)
                                for _ in self.genetic_code.amino_acids], 100)

    def find_best_equimolar_degenerate_codon(self):
        min_cost = 99999
        best_dc = []

        # check all 3315 possible equimolar degenerate codons
        for dc_label in DEGENERATE_CODONS:
            bases = [LiteBase([int(_ in DEGENERATE_NUCLEOTIDE_CODE[d])
                               for _ in 'ACGT'])
                     for d in dc_label]
            dc = LiteCodon(bases, self.genetic_code)
            cost = self.calculate_cost(dc.amino_acid_profile)

            if cost < min_cost:
                min_cost = cost
                best_dc = [dc]
            elif cost == min_cost:
                best_dc.append(dc)

        return best_dc

    def calculate_cost(self, current):
        """
        sum of squared errors
        """
        current = rescale(current, 100)
        return sum([(x - y) ** 2
                    for x, y in zip(current, self.targets)])

    def exhaustive_search(self, precision=20):
        self.sparse_data = []
        bases = []
        for a in range(0, 100+1, precision):
            for c in range(0, 100-a+1, precision):
                for g in range(0, 100-a-c+1, precision):
                    t = 100 - a - c - g
                    bases.append(LiteBase([a, c, g, t]))

        for i in range(len(bases)):
            for j in range(len(bases)):
                for k in range(len(bases)):
                    dc = LiteCodon(bases=[bases[i], bases[j], bases[k]],
                                   genetic_code=self.genetic_code)
                    cost = self.calculate_cost(dc.amino_acid_profile)
                    self.sparse_data.append(
                        dc.bases[0].acgt +
                        dc.bases[1].acgt +
                        dc.bases[2].acgt +
                        [cost]
                        )

        colnames = [
                    'A1', 'C1', 'G1', 'T1',
                    'A2', 'C2', 'G2', 'T2',
                    'A3', 'C3', 'G3', 'T3',
                    'SSE']
        self.sparse_df = pd.DataFrame(self.sparse_data,
                                      columns=colnames)

    def run(self):
        if self.optimization_mode == 'random':
            best_cost = 99999
            for _ in range(self.max_total_samples):
                bases = [LiteBase(get_random_base_profile())
                         for _ in range(3)]
                dc = LiteCodon(bases=bases,
                               genetic_code=self.genetic_code)

                cost = self.calculate_cost(dc.amino_acid_profile)
                if cost < best_cost:
                    print(cost)
                    print(dc.base_profile)
                    best_cost = cost
                    self.best_dc = [dc]
                elif cost == best_cost:
                    self.best_dcs.append(dc)

        for dc in self.best_dcs:
            print(dc.bases[0].acgt, dc.bases[1].acgt, dc.bases[2].acgt)
            print(self.calculate_cost(dc.amino_acid_profile))
            print(dc.amino_acid_profile)


def main():
    opt = NonEquimolarOptimizer()
    opt.run()
    return opt


if __name__ == '__main__':
    opt = main()
