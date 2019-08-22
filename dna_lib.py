VARIETY = 4

NUCLEOTIDES = ['A', 'C', 'G', 'T']


def check_nucleotide(nucleotide):
    return nucleotide in NUCLEOTIDES


def nucleotide_to_number(nucleotide):
    return NUCLEOTIDES.index(nucleotide)


def number_to_nucleotide(number):
    return NUCLEOTIDES[number]


def pattern_to_number(pattern):
    if not pattern:
        return 0

    nucleotide = pattern[-1]

    if not check_nucleotide(nucleotide):
        raise ValueError(format("Invalid nucleotide %s" % nucleotide))

    reminder = nucleotide_to_number(nucleotide)
    quotient = pattern[:-1]

    return pattern_to_number(quotient) * VARIETY + reminder


def number_to_pattern(index, k):
    if k == 1:
        return number_to_nucleotide(index)

    reminder = index % VARIETY
    nucleotide = number_to_nucleotide(reminder)

    prefix_index = index // VARIETY
    prefix_pattern = number_to_pattern(prefix_index, k - 1)
    return prefix_pattern + nucleotide


def format_frequencies_result(arr):
    result = ''
    for item in arr:
        result += format('%d ' % item)
    return result.strip()


def computing_frequencies(text, k):
    n = 4 ** k
    frequency_array = [0] * n
    for i in range(0, len(text) - (k - 1)):
        pattern = text[i: i + k]

        idx = pattern_to_number(pattern)
        frequency_array[idx] += 1

    return frequency_array


def hamming_dist(seq1, seq2):
    dist = abs(len(seq1) - len(seq2))
    if dist:
        return dist

    for s1, s2 in zip(seq1, seq2):
        if s1 != s2:
            dist += 1
    return dist


def immediate_neighbours(pattern):
    neighbours = {pattern}
    for i in range(0, len(pattern)):
        for nucleotide in NUCLEOTIDES:
            new_pattern = pattern[:i] + nucleotide + pattern[i + 1:]
            neighbours.add(new_pattern)
    return neighbours


def neighbors(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return set(NUCLEOTIDES)
    neighborhood = set()
    suffix = pattern[1:]
    suffix_neighbors = neighbors(suffix, d)
    for neighbor in suffix_neighbors:
        if hamming_dist(suffix, neighbor) < d:
            for nucleotide in NUCLEOTIDES:
                neighborhood.add(nucleotide + neighbor)
        else:
            neighborhood.add(pattern[:1] + neighbor)
    return neighborhood


def iterative_neighbors(pattern, d):
    neighborhood = {pattern}
    for j in range(d):
        for neighbor in neighborhood:
            neighborhood = neighborhood.union(immediate_neighbours(neighbor))
    return neighborhood


def exact_neighbors(pattern, d):
    neighborhood = neighbors(pattern, d)
    return {p for p in neighborhood if hamming_dist(pattern, p) == d}


def computing_frequencies_with_mismatches(text, k, d):
    frequency_array = dict()
    frequency_array.setdefault(0)
    for i in range(0, len(text)-k):
        pattern = text[i, k]
        pattern_neighbourhood = neighbors(pattern, d)
        for approximate_pattern in pattern_neighbourhood:
            j = pattern_to_number(approximate_pattern)
            frequency_array[j] += 1
    return frequency_array


def __check_skew_params(genome, i):
    if i < 0:
        raise ValueError("i should be greater or equal to 0")
    if not genome:
        raise ValueError("genome should not be an empty string")


def skew(genome, i):
    """
    Compute Skewi+1(Genome) from Skewi(Genome) according to the nucleotide in position i of Genome. If this
    nucleotide is G, then Skewi+1(Genome) = Skewi(Genome) + 1; if this nucleotide is C, then Skewi+1(Genome)= Skewi(
    Genome) – 1; otherwise, Skewi+1(Genome) = Skewi(Genome).
     :param i: position, calculate skew till
     :param genome: genome to compute the skew characteristics
     :return: result skew
    """
    __check_skew_params(genome, i)
    result = 0
    for nucleotide in genome[:i]:
        if nucleotide == 'G':
            result += 1
        elif nucleotide == 'C':
            result -= 1
    return result
