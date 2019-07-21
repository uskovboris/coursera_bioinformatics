def hamming_dist(seq1, seq2):
    dist = abs(len(seq1) - len(seq2))
    if dist:
        return dist

    for s1, s2 in zip(seq1, seq2):
        if s1 != s2:
            dist += 1
    return dist
