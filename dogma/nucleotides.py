from collections import OrderedDict, defaultdict
from random import choices

from dogma import(rescale,
                  get_frequency_dictionary,

                  STANDARD_NUCLEOTIDES,
                  DEFAULT_NUCLEOTIDE_LABEL,
                  DEGENERATE_NUCLEOTIDE_CODE,
                  DEGENERATE_NUCLEOTIDES,
                  DEGENERATE_NUCLEOTIDE_CODE_COMPOSITION,
                  DEGENERATE_NUCLEOTIDE_CODE_REVERSED,)


class Nucleotide:
    """
    Nucleotides are the building blocks of oligonucleotides.

    Parameters
    ----------
    data: multiple inputs available, see self._process_input()

    Attributes
    ----------
    label: str, len=1
        single letter code for nucleotide
    composition: dictionary
        keys are single letters refering to undegenerate nucleotides
        values are numerical values related to abundance proportion
    members
        sorted list of nonzero nucleotides (as single letter codes)
    proportions
        list of member proportions, aligned with members
    degenerate: Boolean
        whether multiple standard nucleotides are nonzero
    equimolar: Boolean
        whether all nonzero standard nucleotide members have equal proportions

    Usage
    -----
    a = Nucleotide('A')
        a.label --> 'A'
        a.composition --> {'A':1}
        a.members --> ['A']
        a.proportions --> [1]
        a.degenerate --> False
        a.equimolar --> True
    a = Nucleotide({'A':1})


    n = Nucleotide('N')
        a.label --> 'N'
        a.composition --> {'A':1, 'C':1, 'G':1, 'T':1}
        a.members --> ['A', 'C', 'G', 'T']
        a.proportions --> [1, 1, 1, 1]
        a.degenerate --> True
        a.equimolar --> True

    n = Nucleotide({'A':0.1, 'C':0.5, 'G':0.4})
        a.label --> 'v'
        a.composition --> {'A':0.1, 'C':0.5, 'G':0.4}
        a.members --> ['A', 'C', 'G', 'T']
        a.proportions --> [0.1, 0.5, 0.4]
        a.degenerate --> True
        a.equimolar --> False
    """

    def __init__(self, data=None):

        # define self.label and self.composition by processing data input parameter
        self._process_input(data)

        self.members = self.get_members()
        self.proportions = self.get_proportions()

        self.degenerate = self.is_degenerate()
        self.equimolar = self.is_equimolar()


    def _process_input(self, data):
        """
        Various input data formats are supported:
            Nucleotide('A')
            Nucleotide({'A': 3, 'C': 1})
        """

        # ex: Nucleotide('A')
        if isinstance(data, str) and len(data) == 1:

            label = data.upper().replace('U', 'T')
            composition = DEGENERATE_NUCLEOTIDE_CODE_COMPOSITION[label]
        
        # ex: Nucleotide({'A':3, 'C':1})
        elif (isinstance(data, dict) and
              is_valid_nucleotide_string(''.join(data.keys())) and 
              all([isinstance(_, (int, float)) for _ in data.values()])):

            composition = data
            label = nucleotide_composition_to_letter(composition)

        else:
            label = DEFAULT_NUCLEOTIDE_LABEL
            composition = DEGENERATE_NUCLEOTIDE_CODE_COMPOSITION[label]


        self.label = label
        self.composition = composition


    def get_members(self):
        """
        Returns a sorted string of nucleotide members that are have nonzero
        frequency within nucleotide.
        """
        nonzero_members = [k for k, v in self.composition.items() if v > 0]
        return ''.join(sorted(nonzero_members))


    def get_proportions(self):
        """
        Returns a list of proportions aligned with self.members
        """
        return [self.composition[_] for _ in self.members]


    def is_degenerate(self):
        """
        Nucleotide is degenerate if multiple standard nucleotides have nonzero proportions.
        """
        return len(self.members) > 1


    def is_equimolar(self):
        """
        Nucleotide is equimolar if the proportion of each nonzero nucleotide values are equal.
        """
        nonzero_values = [_ for _ in self.composition.values() if _ > 0]
        return max(nonzero_values) == min(nonzero_values)


    def samples(self, k=1):
        """
        Randomly selects `k` members based on Nucleotide composition.
        """
        return [self.sample() for _ in range(k)]


    def sample(self):
        """
        Randomly selects a members based on Nucleotide composition.
        """
        return choices(self.members, self.proportions)[0]


    def copy(self):
        return Nucleotide(self.composition)


    def __str__(self):
        return f'Nucleotide(label={self.label}, composition={self.composition})'


    def __repr__(self):
        return f'Nucleotide(label={self.label}, composition={self.composition})'


def is_nondegenerate_nucleotide_string(n):
    n = n.upper().replace('U', 'T')
    return all([_ in STANDARD_NUCLEOTIDES for _ in n])


def is_degenerate_nucleotide_string(n):
    if not is_nondegenerate_nucleotide_string(n):
        n = n.upper().replace('U', 'T')
        return all([_ in DEGENERATE_NUCLEOTIDES for _ in n])
    return False


def is_valid_nucleotide_string(s):
    """
    Checks if string is composed of valid nucleotide members.

    input is first capitalized and U's are replaced with T's
    each letter must be a standard or degenerate letter
    """
    s = s.upper().replace('U', 'T')
    return all([_ in DEGENERATE_NUCLEOTIDES for _ in s])


def combine_nucleotides(*nucleotides, proportions=None):
    """
    Combines nucleotides, creating a new Nucleotide instance with
    a nucleotide composition as a weighted sum input compositions.

    Usage
    -----
    c = Nucleotide('C')
    t = Nucleotide('T')

    c_plus_t = combine_nucleotides(c, t)
    print(c_plus_t) --> "Nucleotide('N', composition={'A':0, 'C':0.5, 'G':0, 'T':0.5})" 

    c_plus_ttt = combine_nucleotides(c, t, proportions=[1, 3])
    print(c_plus_ttt) --> "Nucleotide('n', composition={'A':0, 'C':0.25, 'G':0, 'T':0.75})" 
    """

    if (proportions is None) or (len(nucleotides) != len(proportions)):
        proportions = [1,] * len(nucleotides)

    data = {}
    for n in STANDARD_NUCLEOTIDES:
        data[n] = 0
        for nuc, p in zip(nucleotides, proportions):
            data[n] += nuc.composition.get(n, 0) * p
    return Nucleotide(rescale(data))


def nucleotide_composition_to_letter(composition):
    """
    Converts dictionary of {nucleotide letter: proportions} pairs to IUPAC degenerate DNA letter.

    Usage:
    c = {'A': 1}
    print(nucleotide_composition_to_letter(c)) --> 'A'

    c = dict(zip('ACGT', [1, 1, 1, 1]))
    print(nucleotide_composition_to_letter(c)) --> 'N'

    c = dict(zip('ACGT', [1, 1, 2, 1]))
    print(nucleotide_composition_to_letter(c)) --> 'n'
    """

    nonzero_nucleotides = ''.join(sorted([n for n, v in composition.items() if v > 0]))
    nonzero_proportions = [composition[n] for n in nonzero_nucleotides]
    equimolar = min(nonzero_proportions) == max(nonzero_proportions)

    letter = DEGENERATE_NUCLEOTIDE_CODE_REVERSED.get(nonzero_nucleotides, DEFAULT_NUCLEOTIDE_LABEL)

    if equimolar:
        return letter
    return letter.lower()
