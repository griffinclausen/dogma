from dogma.extensions.codon_selection import (
    CodonOption,
    CodonSelector
)


def test_codon_selector():
    cs = CodonSelector()
    assert len(cs.table) == 3375
    assert len(cs.filtered_table) == 3375

    cs.must_exclude = ['*']
    cs.filter()
    assert len(cs.filterd_table) = 2351


def main():
    cs = CodonSelector()
    return cs


if __name__ == '__main__':
    main()
