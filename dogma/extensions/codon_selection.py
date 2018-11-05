import pandas as pd
from dogma import (
    Codon,
    GeneticCode,
    DEFAULT_GENETIC_CODE,
    NCBI_CODONS,
    DEGENERATE_CODONS
)


class CodonOption(Codon):
    def __init__(self, data=None, genetic_code=None):
        super().__init__(data, genetic_code)

        self.amino_acid_profile = self.translate()
        self.aa_string = self.get_aa_string()
        self.aas_string = self.get_aas_string()
        self.aa_profile_string = self.get_aa_profile_string()

    def get_aa_string(self):
        aas = sorted([_ for _, v in self.amino_acid_profile.items() if v > 0])
        return ''.join(aas)

    def get_aas_string(self):
        aas = [self.genetic_code[_] for _ in self.members]
        return ''.join(sorted(aas))

    def get_aa_profile_string(self):
        aas = sorted([_ for _, v in self.amino_acid_profile.items() if v > 0])
        output = [f'{aa}:{self.amino_acid_profile[aa]},' for aa in aas]
        return ''.join(output)

    def is_aa_equimolar(self):
        vals = [_ for _ in self.amino_acid_profile.values() if _ > 0]
        return max(vals) == min(vals)


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
        self.amino_acids = genetic_code.amino_acids
        self.codons = genetic_code.codons

        self.must_exclude = must_exclude
        self.must_include = must_include

        self.scope = scope
        self.table = self.make_table()
        self.filter()

    def parse_criteria(self):
        self.excluded_amino_acids = []
        self.excluded_codons = []
        for _ in self.must_exclude:
            if _ in self.amino_acids:
                self.excluded_amino_acids.append(_)
            elif _ in self.codons:
                self.excluded_codons.append(_)

        self.included_amino_acids = []
        self.included_codons = []
        for _ in self.must_include:
            if _ in self.amino_acids:
                self.included_amino_acids.append(_)
            elif _ in self.codons:
                self.included_codons.append(_)

    def make_table(self, dcs=None):
        if dcs is None:
            dcs = [CodonOption(_, self.genetic_code) for _ in self.scope]

        data = {
            'Labels': [_.label for _ in dcs],
            'AAs': [_.get_aas_string() for _ in dcs],
            'uAAs': [_.aa_string for _ in dcs]
            }
        data['nCodons'] = [len(_) for _ in data['AAs']]
        data['nuAAs'] = [len(_) for _ in data['uAAs']]
        data['Stops'] = ['*' in _ for _ in data['uAAs']]
        data['Equimolar'] = [_.is_aa_equimolar() for _ in dcs]

        return pd.DataFrame(data)

    def filter(self):
        self.parse_criteria()

        self.filtered_options = []
        for dc in self.scope:
            c = CodonOption(dc, self.genetic_code)
            data = c.translate()
            PASSES_FILTERS = True

            for _ in self.excluded_codons:
                if _ in c.members:
                    PASSES_FILTERS = False
                    break

            if PASSES_FILTERS:
                for _ in self.excluded_amino_acids:
                    if data.get(_, 0) > 0:
                        PASSES_FILTERS = False
                        break

            if PASSES_FILTERS:
                for _ in self.included_codons:
                    if _ not in c.members:
                        PASSES_FILTERS = False
                        break

            if PASSES_FILTERS:
                for _ in self.included_amino_acids:
                    if data.get(_, 0) == 0:
                        PASSES_FILTERS = False
                        break

            if PASSES_FILTERS:
                self.filtered_options.append(c)

        self.filtered_table = self.make_table(self.filtered_options)


def main():
    CS = CodonSelector(must_exclude=list('*'))
    CS.filter()
    return CS


if __name__ == '__main__':
    cs = main()
