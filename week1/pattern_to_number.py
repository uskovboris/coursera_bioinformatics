VARIETY=4

ADENINE='A'
TIMIN='T'
GUANIN='G'
CYTOSINE='C'

def check_nucleotide(nucleotide):
    return nucleotide in [ADENINE, TIMIN, GUANIN, CYTOSINE]

def nucleotide_to_number(nucleotide):
    if nucleotide == ADENINE:
        return 1
    elif nucleotide == TIMIN:
        return 2
    elif nucleotide == GUANIN:
        return 3
    elif nucleotide == CYTOSINE:
        return 4

def pattern_to_number(pattern):
    nucleotide = pattern[-2:-1]
    print('nucleotide: %s' % nucleotide)

    if not check_nucleotide(nucleotide):
        raise ValueError(format("Invalid nucliotide %s" % nucleotide))

    reminder = nucleotide_to_number(nucleotide)
    qantity = pattern[:-1]

    if len(qantity)>1:
        return pattern_to_number(qantity) * VARIETY + reminder
    else:
        return nucleotide_to_number(qantity)


print(pattern_to_number("AAA"))
print(pattern_to_number("AAB"))
print(pattern_to_number("AAT"))
print(pattern_to_number("AAG"))