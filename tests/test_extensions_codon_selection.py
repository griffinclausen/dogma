from dogma.extensions import (
    CodonSelector
)


def test_codon_selector():
    cs = CodonSelector()
    assert len(cs.table) == 3375
    assert len(cs.filtered_table) == 3375

    cs.must_exclude = ['*']
    cs.filter()
    assert len(cs.filtered_table) == 2351
