from dogma.utils import (rescale,
                         get_frequency_dictionary)


def test_rescale_with_list_input():

    d = [1, 1, 1, 1]
    assert isinstance(rescale(d), list)
    assert len(rescale(d)) == len(d)
    assert sum(rescale(d)) == 1
    assert sum(rescale(d, 1)) == 1
    assert sum(rescale(d, 2)) == 2
    assert min(rescale(d, 17)) == max(rescale(d, 17))
    assert rescale(d, 4) == d

    d = [1.0, 21, 0, 5]
    assert isinstance(rescale(d), list)
    assert len(rescale(d)) == len(d)
    assert sum(rescale(d)) == 1
    assert sum(rescale(d, 1)) == 1
    assert sum(rescale(d, 2)) == 2
    assert rescale(d, 27) == d

    d = [2]
    assert isinstance(rescale(d), list)
    assert len(rescale(d)) == len(d)
    assert sum(rescale(d)) == 1
    assert sum(rescale(d, 1)) == 1
    assert sum(rescale(d, 2)) == 2
    assert rescale(d, 2) == d


def test_rescale_with_dict_input():

    d = dict(zip('abcd', [1, 1, 1, 1]))
    assert isinstance(rescale(d), dict)
    assert len(rescale(d)) == len(d)
    assert sum(rescale(d).values()) == 1
    assert sum(rescale(d, 1).values()) == 1
    assert sum(rescale(d, 2).values()) == 2
    assert min(rescale(d, 17).values()) == max(rescale(d, 17).values())
    assert rescale(d, 4) == d

    d = dict(zip('abcd', [1.0, 21, 0, 5]))
    assert isinstance(rescale(d), dict)
    assert len(rescale(d)) == len(d)
    assert sum(rescale(d).values()) == 1
    assert sum(rescale(d, 1).values()) == 1
    assert sum(rescale(d, 2).values()) == 2
    assert rescale(d, 27) == d

    d = {'a': 2}
    assert isinstance(rescale(d), dict)
    assert len(rescale(d)) == len(d)
    assert sum(rescale(d).values()) == 1
    assert sum(rescale(d, 1).values()) == 1
    assert sum(rescale(d, 2).values()) == 2
    assert rescale(d, 2) == d


def test_get_frequency_dictionary_with_list_input():

    x = list('aaabbc')
    y = get_frequency_dictionary(x)
    assert y == {'a': 3, 'b': 2, 'c': 1}

    x = ['aaa', 'bbb', 'aaa', 'bbb', 'aaa', 'ccc']
    y = get_frequency_dictionary(x)
    assert y == {'aaa': 3, 'bbb': 2, 'ccc': 1}


def test_get_frequency_dictionary_with_string_input():

    x = 'aaabbc'
    y = get_frequency_dictionary(x)
    assert y == {'a': 3, 'b': 2, 'c': 1}

    x = ['aaa', 'bbb', 'aaa', 'bbb', 'aaa', 'ccc']
    y = get_frequency_dictionary(x)
    assert y == {'aaa': 3, 'bbb': 2, 'ccc': 1}
