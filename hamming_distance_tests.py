import unittest
from hamming_dist import *


class ReadFileTest(unittest.TestCase):

    def test_hamming_dist_0(self):
        """ Hamming distance test """
        self.assertEqual(hamming_dist([1, 2, 3], [1, 2, 3]), 0)

    def test_hamming_dist_1(self):
        self.assertEqual(hamming_dist([1, 2, 3], [1, 4, 3]), 1)
        self.assertEqual(hamming_dist([1, 2, 3], [1, 4, 5]), 2)
        self.assertEqual(hamming_dist([1, 2, 3], [2, 4, 5]), 3)
        self.assertEqual(hamming_dist([1, 4, 3], [1, 2, 3]), 1)
        self.assertEqual(hamming_dist([1, 4, 5], [1, 2, 3]), 2)
        self.assertEqual(hamming_dist([1, 4], [1, 2, 3]), 2) # ???
        self.assertEqual(hamming_dist([1], [1, 2, 3]), 2) # ???

    def test_hamming_dist_strings(self):
        self.assertEqual(hamming_dist("GGGCCGTTGGT", "GGACCGTTGAC"), 3)


if __name__ == "__main__":
    unittest.main()