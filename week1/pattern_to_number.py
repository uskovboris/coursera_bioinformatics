VARIETY = 4


def check_nucleotide(nucleotide):
    return nucleotide in ['A', 'C', 'G', 'T']


def nucleotide_to_number(nucleotide):
    if nucleotide == 'A':
        return 0
    elif nucleotide == 'C':
        return 1
    elif nucleotide == 'G':
        return 2
    elif nucleotide == 'T':
        return 3


def pattern_to_number(pattern):
    if not pattern:
        return 0

    nucleotide = pattern[-1]

    if not check_nucleotide(nucleotide):
        raise ValueError(format("Invalid nucleotide %s" % nucleotide))

    reminder = nucleotide_to_number(nucleotide)
    quantity = pattern[:-1]

    return pattern_to_number(quantity) * VARIETY + reminder


print(pattern_to_number("TCCTAGTAGTCAAGCTC"))
