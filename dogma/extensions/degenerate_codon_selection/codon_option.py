from dogma import Codon


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
