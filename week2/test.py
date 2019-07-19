from hamming_dist import *


def main():
    print(hamming_dist([1, 2, 3], [1, 2, 3]))

    print(hamming_dist([1, 2, 3], [1, 4, 3]))
    print(hamming_dist([1, 2, 3], [1, 4, 5]))
    print(hamming_dist([1, 2, 3], [2, 4, 5]))

    print(hamming_dist([1, 4, 3], [1, 2, 3]))
    print(hamming_dist([1, 4, 5], [1, 2, 3]))

    print(hamming_dist([1, 4], [1, 2, 3]))
    print(hamming_dist([1], [1, 2, 3]))


if __name__ == '__main__':
    main()
