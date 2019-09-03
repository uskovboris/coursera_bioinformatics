from dna_lib import *
from fs_helpers import *

def hamming_dist_in(cur_dna_chunk, pattern, k, d):
    for j in range(0, len(cur_dna_chunk) - k + 1):
        cur_pattern = cur_dna_chunk[j:j + k]
        if hamming_dist(pattern, cur_pattern) <= d:
            return True
    else:
        return False


def motifs_enumeration(dna, k, d):
    motifs = []
    complete_dna = "".join(dna)
    for i in range(0, len(complete_dna) - k + 1):
        pattern = complete_dna[i:i + k]
        neighborhood = neighbors(pattern, d)
        for approximate_pattern in neighborhood:
            if all(hamming_dist_in(l, approximate_pattern, k, d) for l in dna):
                if approximate_pattern not in motifs:
                    motifs.append(approximate_pattern)
    return sorted(motifs)


if __name__ == '__main__':
    lines = read_lines(r"C:\Users\Boris\Downloads\dataset_156_8.txt")
    k, d = (int(p) for p in lines[0].split())
    dna = lines[1:]

    motifs = motifs_enumeration(dna, k, d)
    print(' '.join(motifs))
