def count_pattern(chromosome, pattern):
    pattern_len = len(pattern)
    count = 0
    for i in range(len(chromosome)):
        if chromosome[i:i + pattern_len] == pattern:
            count += 1
    return count


def frequent_word(chromosome, pattern_len):
    counts_dict = dict()
    max_pattern_count = 0
    for i in range(len(chromosome) - pattern_len):
        current_pattern = chromosome[i:i + pattern_len]
        counts_dict[current_pattern] = count_pattern(chromosome, current_pattern, pattern_len)
        print(current_pattern, counts_dict[current_pattern])
        if max_pattern_count < counts_dict[current_pattern]:
            max_pattern_count = counts_dict[current_pattern]

    return list(pattern for pattern, cnt in counts_dict.items() if cnt == max_pattern_count)


def pattern_pos(chromosome, pattern):
    result = ""
    pattern_len = len(pattern)
    for i in range(len(chromosome) - pattern_len):
        if chromosome[i:i + pattern_len] == pattern:
            result += str(i) + ' '
    return result
