import sys
import collections
import numpy

VARIETY = 4

DNA_NUCLEOTIDES = ['A', 'C', 'G', 'T']
RNA_NUCLEOTIDES = ['A', 'C', 'G', 'U']


def check_nucleotide(nucleotide, RNA=False):
    return nucleotide in DNA_NUCLEOTIDES if not RNA else RNA_NUCLEOTIDES


def nucleotide_to_number(nucleotide, RNA=False):
    if not RNA:
        if check_nucleotide(nucleotide):
            return DNA_NUCLEOTIDES.index(nucleotide)
        else:
            raise ValueError("Invalid for DNA nucleotide '{}'".format(nucleotide))

    elif RNA:
        if check_nucleotide(nucleotide, RNA=True):
            return RNA_NUCLEOTIDES.index(nucleotide)
        else:
            raise ValueError("Invalid for DNA nucleotide '{}'".format(nucleotide))


def number_to_nucleotide(number, RNA=False):
    return DNA_NUCLEOTIDES[number] if not RNA else RNA_NUCLEOTIDES[number]


def pattern_to_number(pattern):
    if not pattern:
        return 0

    nucleotide = pattern[-1]

    if not check_nucleotide(nucleotide):
        raise ValueError("Invalid nucleotide {}".format(nucleotide))

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


def format_list_result(arr):
    result = ''
    for item in arr:
        result += '{} '.format(item)
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


def hamming_dist_in(cur_dna_chunk, pattern, k, d):
    for j in range(0, len(cur_dna_chunk) - k + 1):
        cur_pattern = cur_dna_chunk[j:j + k]
        if hamming_dist(pattern, cur_pattern) <= d:
            return True
    else:
        return False


def immediate_neighbours(pattern):
    neighbours = {pattern}
    for i in range(0, len(pattern)):
        for nucleotide in DNA_NUCLEOTIDES:
            new_pattern = pattern[:i] + nucleotide + pattern[i + 1:]
            neighbours.add(new_pattern)
    return list(neighbours)


def neighbors(pattern, d):
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return DNA_NUCLEOTIDES
    neighborhood = set()
    suffix = pattern[1:]
    suffix_neighbors = neighbors(suffix, d)
    for neighbor in suffix_neighbors:
        if hamming_dist(suffix, neighbor) < d:
            for nucleotide in DNA_NUCLEOTIDES:
                neighborhood.add(nucleotide + neighbor)
        else:
            neighborhood.add(pattern[:1] + neighbor)
    return list(neighborhood)


def iterative_neighbors(pattern, d):
    neighborhood = {pattern}
    for j in range(d):
        for neighbor in neighborhood:
            neighborhood = neighborhood.union(immediate_neighbours(neighbor))
    return list(neighborhood)


def exact_neighbors(pattern, d):
    neighborhood = neighbors(pattern, d)
    return [p for p in neighborhood if hamming_dist(pattern, p) == d]


def computing_frequencies_with_mismatches(text, k, d, with_complements=False):
    n = 4 ** k
    frequency_array = [0] * n
    for i in range(0, len(text)-k+1):
        pattern = text[i: i+k]
        pattern_neighbourhood = neighbors(pattern, d)
        if with_complements:
            pattern_reverse_complement = reverse_complement(pattern)
            pattern_neighbourhood.extend(neighbors(pattern_reverse_complement, d))
        for approximate_pattern in pattern_neighbourhood:
            pattern_hash = pattern_to_number(approximate_pattern)
            frequency_array[pattern_hash] += 1
    return frequency_array


def reverse_complement(template_strand, RNA=False):
    complement_strand = ""
    for nucleotide in template_strand:
        if nucleotide == 'A':
            complement_strand = 'T' + complement_strand
        elif nucleotide == 'T':
            if RNA:
                complement_strand = 'U' + complement_strand
            else:
                complement_strand = 'A' + complement_strand
        elif nucleotide == 'G':
            complement_strand = 'C' + complement_strand
        elif nucleotide == 'C':
            complement_strand = 'G' + complement_strand
        elif nucleotide == 'U':
            if RNA:
                complement_strand = 'T' + complement_strand
            else:
                raise ValueError("Invalid for DNA nucleotide '{}'".format(nucleotide))
    return complement_strand


def computing_frequent_patterns_with_mismatches(text, k, d, with_complements=False):
    frequency_dict = computing_frequencies_with_mismatches(text, k, d, with_complements)

    max_frequency = max(frequency_dict)
    return [number_to_pattern(i, k) for i, val in enumerate(frequency_dict) if (val == max_frequency)]


def approximate_patterns_count(text, pattern, d):
    count = 0
    pattern_len = len(pattern)
    for i in range(len(text)-pattern_len+1):
        current_pattern = text[i:i+pattern_len]
        if hamming_dist(current_pattern, pattern) <= d:
            count += 1
    return count


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


def __check_find_skew_minimums_params(genome):
    if not genome:
        raise ValueError("genome should not be an empty string")


def find_skew_minimums(genome):
    __check_find_skew_minimums_params(genome)
    skews = dict({0: 0})
    min_skew = len(genome)
    for i, nucleotide in enumerate(genome):
        if nucleotide == 'G':
            skews[i + 1] = (skews[i] + 1)
        elif nucleotide == 'C':
            skews[i + 1] = (skews[i] - 1)
        else:
            skews[i + 1] = skews[i]

        if min_skew > skews[i + 1]:
            min_skew = skews[i + 1]

    return [skew_i for skew_i in range(len(skews)) if skews[skew_i] == min_skew], skews


def __pattern_positions_internal(pattern, chromosome, matcher):
    result = []
    pattern_len = len(pattern)
    for i in range(len(chromosome) - pattern_len + 1):
        if matcher(chromosome[i:i + pattern_len], pattern):
            result.append(i)
    return result


def pattern_positions(pattern, chromosome):
    return __pattern_positions_internal(pattern, chromosome, lambda str1, str2: str1 == str2)


def approximate_pattern_positions(pattern, chromosome, d):
    return __pattern_positions_internal(pattern, chromosome, lambda str1, str2: hamming_dist(str1, str2) <= d)


def entropy(vector):
    """
    A measure of the uncertainty of a probability distribution (p1, …, pN)
    :param vector: vector to calculate an entropy
    """
    return -sum([p * numpy.math.log(p, 2) for p in vector if p > 0])


def matrix_entropy(motif_matrix):
    """Calculates entropy of motifs matrix as a sum of entropies of the matrix columns"""
    result = 0
    motif_matrix = numpy.array(motif_matrix)
    columns = motif_matrix.T
    for column in columns:
        counts = collections.Counter(column)
        frequencies = [cnt/len(column) for cnt in counts.values()]
        result += entropy(frequencies)
    return result


def motifs_enumeration(dna, k, d):
    """
    Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif if it appears in every string from Dna
    with at most d mismatches. For example, the implanted 15-mer in the strings above represents a (15,4)-motif.
    :param dna: the DNA strand to find motifs
    :param k: expected motif length
    :param d: hamming distance upper limit
    :return: the list of found motifs
    """
    motifs = []
    complete_dna = "".join(dna)
    for i in range(0, len(complete_dna) - k + 1):
        pattern = complete_dna[i:i + k]
        neighborhood = neighbors(pattern, d)
        for approximate_pattern in neighborhood:
            if all(hamming_dist_in(cur_dna_chunk, approximate_pattern, k, d) for cur_dna_chunk in dna):
                if approximate_pattern not in motifs:
                    motifs.append(approximate_pattern)
    return sorted(motifs)


def distance_between_pattern_and_strings(pattern, dna):
    k = len(pattern)
    distance = 0
    for cur_dna_chunk in dna:
        min_hamming_dist = 100000
        for i in range(0, len(cur_dna_chunk) - k + 1):
            cur_pattern = cur_dna_chunk[i:i + k]
            cur_pattern_hamming_dist = hamming_dist(pattern, cur_pattern)
            if cur_pattern_hamming_dist < min_hamming_dist:
                min_hamming_dist = cur_pattern_hamming_dist
        distance += min_hamming_dist
    return distance


def median_string(dna, k):
    n = 4 ** k
    patterns = [number_to_pattern(i, k) for i in range(0, n)]

    median = ""
    distance = 100000
    for pattern in patterns:
        cur_pattern_distance = distance_between_pattern_and_strings(pattern, dna)
        if distance > cur_pattern_distance:
            distance = cur_pattern_distance
            median = pattern
    return median
