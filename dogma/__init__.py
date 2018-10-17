"Collection of objects related by the central dogma of biology."
__version__ = '0.0.1'
name = 'dogma'


from dogma.nbci import (
    NBCI_GENETIC_CODE_NAMES,
    NCBI_NUCLEOTIDES,
    NCBI_CODONS,
    NCBI_AMINO_ACIDS,
)

from dogma.utils import (
    DEFAULT_DECIMAL_PRECISION,
    STANDARD_NUCLEOTIDES,
    DEGENERATE_NUCLEOTIDE_CODE,
    DEGENERATE_NUCLEOTIDES,
    DEGENERATE_NUCLEOTIDE_CODE_REVERSED,
    DEGENERATE_NUCLEOTIDE_CODE_COMPOSITION,
    NUCLEOTIDE_BASE_PAIRS,
    DEGENERATE_NUCLEOTIDE_PAIRS,
    DEFAULT_NUCLEOTIDE_LABEL,
    DEFAULT_OLIGONUCLEOTIDE_LABEL,
    STANDARD_CODONS,
    DEGENERATE_CODONS,
    DEFAULT_CODON_LABEL,
    STANDARD_AMINO_ACIDS,
    DEFAULT_AMINO_ACID_LABEL,
    STOP_LABEL,
    DEFAULT_RESIDUE_LABEL,
    rescale,
    get_frequency_dictionary,
    get_random_oligonucleotide
)

from dogma.genetic_codes import (
    GeneticCode,
    DEFAULT_GENETIC_CODE,
    translate)

from dogma.nucleotides import (
    Nucleotide,
    combine_nucleotides,
    is_valid_nucleotide_string
)

from dogma.codons import (
    Codon,
    degenerate_codon_string_to_standard_members,
    nucleotide_string_to_codons,
    combine_codons,
    combine_codon_labels,
    merge_codons
)

from dogma.oligonucleotides import (
    Oligonucleotide,
    reverse_complement
)

from dogma.amino_acids import AminoAcid

from dogma.proteins import (
    Protein,
    calculate_protein_degeneracy
)
