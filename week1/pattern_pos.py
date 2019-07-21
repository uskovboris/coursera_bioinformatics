import sys


def show_usage():
    print("Usage:")
    print("python pattern_count.py <dataset_file>")


def pattern_pos(pattern, chromosome):
    result = ""
    pattern_len = len(pattern)
    for i in range(len(chromosome) - pattern_len + 1):
        if chromosome[i:i + pattern_len] == pattern:
            if result != "":
                result += " "
            result += str(i)
    return result


assert("1 3 9" == pattern_pos("ATAT", "GATATATGCATATACTT"))
assert("4" == pattern_pos("ACAC", "TTTTACACTTTTTTGTGTAAAAA"))
assert("0 46 51 74" == pattern_pos("AAA", "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTACACAACATCCAT"))
assert("88 92 98 132" == pattern_pos("TTT", "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT"))
assert("0 2 4" == pattern_pos("ATA", "ATATATA"))


def main():
    if len(sys.argv) != 2:
        show_usage()

    dataset_file = sys.argv[1]

    with open(dataset_file) as f:
        dataset = f.readlines()

    if len(dataset) != 2:
        print('Dataset should contains 2 lines')
        sys.exit(1)

    chromosome = dataset[1].strip()
    pattern = dataset[0].strip()

    print("chromosome: %s" % chromosome)
    print("pattern: %s" % pattern)

    output = pattern_pos(chromosome, pattern)
    print(output)
    with open('./output.txt', "w") as f:
        f.write(output)

#
# if __name__ == '__main__':
#     main()
