import sys


def show_usage():
    print("Usage:")
    print("python pattern_count.py <pattern> <chromosome_file>")


def pattern_pos(chromosome, pattern):
    result = ""
    pattern_len = len(pattern)
    for i in range(len(chromosome) - pattern_len):
        if chromosome[i:i + pattern_len] == pattern:
            # print("%d: '%s'" % (i , chromosome[i:i + pattern_len]))
            result += str(i) + ' '
    return result


def main():
    if len(sys.argv) != 3:
        show_usage()

    chromosome_file = sys.argv[2]

    with open(chromosome_file) as f:
        chromosome = f.read()

    pattern = sys.argv[1].strip()

    print("chromosome: %s" % chromosome)
    print("pattern: %s" % pattern)

    output = pattern_pos(chromosome, pattern)
    print(output)
    with open('./output.txt', "w") as f:
        f.write(output)


if __name__ == '__main__':
    main()
