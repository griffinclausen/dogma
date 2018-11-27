import pandas as pd
from dogma import (
    DEFAULT_GENETIC_CODE,
    DEGENERATE_CODONS
)

from .codon_option import CodonOption


class CodonSelector:
    """
    Container for the selection of degenerate codons within a genetic code that
    fulfill filtering criteria, such as exlcuding stop codons, rare codons, or
    must include specific amino acids.
    """

    def __init__(self,
                 genetic_code=DEFAULT_GENETIC_CODE,
                 must_include=[],
                 must_exclude=[],
                 scope=DEGENERATE_CODONS):

        self.genetic_code = genetic_code
        self.amino_acids = list(genetic_code.amino_acids)
        self.codons = genetic_code.codons

        self.must_exclude = must_exclude
        self.must_include = must_include

        self.scope = scope
        self.table = self._make_table()
        self.filter()

    def _parse_criteria(self):
        self.excluded_amino_acids = [_ for _ in self.must_exclude
                                     if _ in self.amino_acids]
        self.excluded_codons = [_ for _ in self.must_exclude
                                if _ in self.codons]

        self.included_amino_acids = [_ for _ in self.must_include
                                     if _ in self.amino_acids]
        self.included_codons = [_ for _ in self.must_include
                                if _ in self.codons]

        # self.excluded_amino_acids = []
        # self.excluded_codons = []
        # for _ in self.must_exclude:
        #     if _ in self.amino_acids:
        #         self.excluded_amino_acids.append(_)
        #     elif _ in self.codons:
        #         self.excluded_codons.append(_)

        # self.included_amino_acids = []
        # self.included_codons = []
        # for _ in self.must_include:
        #     if _ in self.amino_acids:
        #         self.included_amino_acids.append(_)
        #     elif _ in self.codons:
        #         self.included_codons.append(_)

    def _make_table(self, dcs=None):
        if dcs is None:
            dcs = [CodonOption(_, self.genetic_code) for _ in self.scope]

        data = {
            'Labels': [_.label for _ in dcs],
            'AAs': [_.get_aas_string() for _ in dcs],
            'uAAs': [_.aa_string for _ in dcs]
            }
        data['Num_Codons'] = [len(_) for _ in data['AAs']]
        data['Num_uAAs'] = [len(_) for _ in data['uAAs']]
        data['Stops'] = ['*' in _ for _ in data['uAAs']]
        data['EquimolarAA'] = [_.is_aa_equimolar() for _ in dcs]

        return pd.DataFrame(data)

    def filter(self):
        self._parse_criteria()

        self.filtered_options = []
        for dc in self.scope:
            c = CodonOption(dc, self.genetic_code)
            data = c.translate()
            PASSES_FILTERS = True

            # exclude option if any excluded codon encountered
            for _ in self.excluded_codons:
                if _ in c.members:
                    PASSES_FILTERS = False
                    break

            # exclude option if any excluded amino acid encountered
            if PASSES_FILTERS:
                for _ in self.excluded_amino_acids:
                    if data.get(_, 0) > 0:
                        PASSES_FILTERS = False
                        break

            # exclude option if any required codons not encountered
            if PASSES_FILTERS:
                for _ in self.included_codons:
                    if _ not in c.members:
                        PASSES_FILTERS = False
                        break

            # exclude option if any required amino acids not encountered
            if PASSES_FILTERS:
                for _ in self.included_amino_acids:
                    if data.get(_, 0) == 0:
                        PASSES_FILTERS = False
                        break

            if PASSES_FILTERS:
                self.filtered_options.append(c)

        self.filtered_table = self._make_table(self.filtered_options)


def main():
    CS = CodonSelector(must_exclude=list('*'))
    CS.filter()
    return CS


if __name__ == '__main__':
    cs = main()
