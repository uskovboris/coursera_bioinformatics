NUCLEOTIDES = ['A', 'C', 'G', 'T']


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
            new_pattern = pattern[:i] + nucleotide + pattern[i+1:]
            neighbours.add(new_pattern)
    return neighbours


# print(immediate_neighbours('ACGT'))


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
