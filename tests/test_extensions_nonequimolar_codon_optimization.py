from dogma.extensions import (
    NonEquimolarOptimizer,
    get_random_base_profile,
    get_random_codon_profile,
    calculate_num_base_combinations,
    calculate_num_codon_combinations
)


def test_nonequimolar_codon_optimizer():
    opt = NonEquimolarOptimizer()
    return opt


def test_calculate_num_codon_combinations():
    for p in [100, 50, 25, 20, 10, 5, 4, 2, 1]:
        print(p, calculate_num_codon_combinations(p))


def test_calculate_num_base_combinations():
    for p in [100, 50, 25, 20, 10, 5, 4, 2, 1]:
        print(p, calculate_num_base_combinations(p))
    print('0.1', calculate_num_base_combinations(1, 1000))


def test_get_random_base_profile():
    for _ in range(10):
        b = get_random_base_profile()
        print(b)


def test_get_random_codon_profile():
    for _ in range(10):
        c = get_random_codon_profile()
        print(c)


def tests():
    # print('testing random base profiles')
    test_get_random_base_profile()

    # print('testing random codon profiles')
    test_get_random_codon_profile()

    # print('testing base combinations')
    test_calculate_num_base_combinations()

    # print('testing codon combinations')
    test_calculate_num_codon_combinations()


def main():
    tests()


if __name__ == '__main__':
    main()
