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


pattern = "TCCTAGTAGTCAAGCTC"
number = pattern_to_number(pattern)
restored_pattern = number_to_pattern(number, len(pattern))

print(number)
print(restored_pattern)

assert (restored_pattern == pattern)
