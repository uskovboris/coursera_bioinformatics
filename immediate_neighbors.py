NUCLEOTIDES = ['A', 'C', 'G', 'T']


def immediate_neighbours(pattern):
    neighbours = {pattern}
    for i in range(0, len(pattern)):
        for nucleotide in NUCLEOTIDES:
            new_pattern = pattern[:i] + nucleotide + pattern[i+1:]
            neighbours.add(new_pattern)
    return neighbours


print(immediate_neighbours('ACGT'))
