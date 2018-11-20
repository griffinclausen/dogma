import os

from dogma.extensions.spiked_oligonucleotides.mixed_oligonucleotides import (
    SpikedOligo)


here = os.path.abspath(os.path.dirname(__file__))


class OligoMixer:
    """
    Container for generating SpikedOligo objects based on from a 'template'
    Oligonucleotide and a spiked-in 'added_oligo'.
    """

    def __init__(self, template, added_oligo):
        self.template = template
        self.added_oligo = added_oligo

    def sweep(self, start=0, stop=101, step=2, scale=100):
        """
        Creates a list of mixed Oligonucleotides generated from the combination
        of a range of proportions (default ratios 0 thtough 100% by 1%)
        """
        self.ratios = [_ / scale for _ in range(start, stop, step)]
        self.make_spiked_oligos(self.ratios)

    def make_spiked_oligo(self, ratio):
        assert ratio >= 0
        assert ratio <= 1
        return SpikedOligo(self.template, self.added_oligo, ratio)

    def make_spiked_oligos(self, ratios):
        self.ratios = ratios
        self.spiked_oligos = [self.make_spiked_oligo(ratio)
                              for ratio in ratios]
        return self.spiked_oligos
