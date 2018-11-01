# dogma

<b>d</b>egenerate<br>
<b>o</b>ligonucleotide<br>
<b>g</b>eneration and<br>
<b>m</b>athematical<br>
<b>a</b>ssessment<br>


# Installation
## Pipenv (or pip)
	>> pipenv install dogma
    >> pip install dogma

## or Github
    >> git clone https://github.com/griffinclausen/dogma.git
    >> cd dogma
    >> pipenv install .

# Basic Usage
## Interactive Python session
    >> import dogma
    >> from dogma import GeneticCode, Oligonucleotide
    >> supE = GeneticCode(1, {'TAG': 'Q'})
    >> oligo = Oligonucleotide('NNK'*7, supE)
    >> samples = oligo.samples(100)
    >> print(samples)

# Contact Information
dogma_maintainer@gmail.com
