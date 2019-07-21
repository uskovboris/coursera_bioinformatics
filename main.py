import sys
from dna_lib import *


def show_usage():
    print("Usage:")
    print("python pattern_count.py <dataset_file>")


def main():
    if len(sys.argv) != 2:
        show_usage()

    dataset_file = sys.argv[1]

    with open(dataset_file) as f:
        dataset = f.readlines()

    if len(dataset) != 2:
        print('Dataset should contains 2 lines')
        sys.exit(1)

    chromosome = dataset[0]
    word_len = int(dataset[1])
    frequent_words = frequent_word(chromosome, word_len)

    print("chromosome: %s" % chromosome)
    print("frequent words:")
    for word in frequent_words:
        print(word)


if __name__ == '__main__':
    main()
