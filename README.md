# dogma

degenerate<br>
 oligonucleotide<br>
  generation and<br>
   mathematical<br>
    assessment<br>


# Installation
## PyPi
    >> pip install dogma

## Github
    >> git clone
    >> cd
    >> setup.py .dogma


# Basic Usage
## Interactive Python session
    >> import dogma
    >> from dogma import GeneticCode, Oligonucleotide
    >> supE = GeneticCode(1, {'TAG': 'Q'})
    >> oligo = Oligonucleotide('NNK'*7, supE)
    >> samples = oligo.samples(100)
    >> print(samples)

## Launch GUI
    >> dogma_gui

## Launch dogma web application on localhost
    >> dogma_app


# Contact Information
Griffin Clausen
