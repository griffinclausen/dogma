# dogma

Collection of object-oriented entities related by the central dogma of biology.

**d**egenerate<br>
**o**ligonucleotide<br>
**g**eneration and<br>
**m**athematical<br>
**a**ssessment<br>


# Installation
## using Pipenv
    >> pipenv install dogma
## or using pip
    >> pip install dogma
## or Github
    >> git clone https://github.com/griffinclausen/dogma.git
    >> cd dogma
    >> pipenv install .

# Basic Usage
    >> from dogma import GeneticCode, Oligonucleotide
    >> supE = GeneticCode(1, {'TAG': 'Q'})
    >> oligo = Oligonucleotide('NNK', supE)
    >> samples = oligo.samples(100)
    >> print(samples)

# Contact Information
dogma_maintainer@gmail.com
