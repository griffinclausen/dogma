from random import shuffle, randint

from dogma import DEGENERATE_NUCLEOTIDE_CODE


def get_random_base_profile(total=100, num_bases=4):
    data = []
    for _ in range(num_bases - 1):
        current_total = sum(data)
        data.append(randint(0, total-current_total))
    current_total = sum(data)
    data.append(total - current_total)
    shuffle(data)
    return data


def get_random_base_profile_within_scope(total=100, scope='N'):
    """
    Returns a list of random integers, of length 'num_bases' and sum of 'total'
    If num_bases is 4, the assumed order is 'ACGT' and the random
    """
    a, c, g, t = 0, 0, 0, 0
    nonzero_bases = DEGENERATE_NUCLEOTIDE_CODE[scope]
    num_bases = len(nonzero_bases)
    random_bases = get_random_base_profile(total, num_bases)

    if 'A' in DEGENERATE_NUCLEOTIDE_CODE[scope]:
        a = random_bases.pop()
    if 'C' in DEGENERATE_NUCLEOTIDE_CODE[scope]:
        c = random_bases.pop()
    if 'G' in DEGENERATE_NUCLEOTIDE_CODE[scope]:
        g = random_bases.pop()
    if 'T' in DEGENERATE_NUCLEOTIDE_CODE[scope]:
        t = random_bases.pop()
    return [a, c, g, t]


def get_random_codon_profile(totals=100):
    return [get_random_base_profile(totals) for _ in range(3)]


def calculate_num_base_combinations(precision=1,
                                    limit=100,
                                    num_bases=4,
                                    spent=0):
    """
    Returns the number of unique base compositions (A, C, G, T) for a given
    DNA synthesis precision. A precisions of 1 corresponds to 1% control over
    base composition (ex: A, C, G, T = 99, 1, 0, 0)

    For 4 base, N = sum[i*(i+1)/2] for i = 1 to n+1, where n = 100 / precision

    This function implements a recurvise solution.
    """
    assert isinstance(precision, int)

    if num_bases > 2:
        return int(sum(
            [calculate_num_base_combinations(precision,
                                             limit,
                                             num_bases-1,
                                             spent+i)
             for i in range(0, limit-spent+1, precision)]
            ))

    elif num_bases == 2:
        return int((limit - spent) / precision + 1)

    elif num_bases == 1:
        return 1


def calculate_num_codon_combinations(precision=1):
    """
    Returns the number of unique codon compositions (three bases comprised of
    A, C, G, T each with base discrete values between 0 and 100 dependant on
    the DNA synthesis precision)
    """
    return calculate_num_base_combinations(precision) ** 3
