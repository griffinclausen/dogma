from dogma import (
    STANDARD_AMINO_ACIDS,
    GeneticCode
)

from dogma.utils import (
    rescale
)

from dogma.extensions.nonequimolar_codon_optimization.lite_dogma import (
    LiteGeneticCode
)

supE = GeneticCode(1, {'TAG': 'Q'})
DEFAULT_GENETIC_CODE = LiteGeneticCode(supE)
DEFAULT_AMINO_ACID_TARGETS = rescale(
    {aa: 1 for aa in STANDARD_AMINO_ACIDS}, 100)
DEFAULT_OPTIMIZATION_MODE = 'random'
DEFAULT_MAX_TOTAL_SAMPLES = 100_000

# Grouping default parameters into dictionary object
DEFAULTS = {
    'genetic_code': DEFAULT_GENETIC_CODE,
    'amino_acid_targets': DEFAULT_AMINO_ACID_TARGETS,
    'optimization_mode': DEFAULT_OPTIMIZATION_MODE,
    'max_total_samples': DEFAULT_MAX_TOTAL_SAMPLES
}
