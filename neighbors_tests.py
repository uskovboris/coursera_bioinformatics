import unittest
from neighbors import *


class NeighboursTest(unittest.TestCase):

    def test_neighbours_d_equals_0(self):
        """d=0 test"""
        self.assertSetEqual(neighbors('ACG', 0), {"ACG"})

    def test_neighbours_len_equals_1(self):
        """1 nucleotide len DNA test"""
        self.assertSetEqual(neighbors('A', 1), {"A", "T", "C", "G"})

    def test_neighbours(self):
        """ neighbours test """
        self.assertSetEqual(neighbors('ACG', 1), {'CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 'ACA', 'ACC', 'ACT', 'ACG'})
        self.assertSetEqual(neighbors('CAA', 1), {'CAA', 'CCA', 'CGA', 'CTA', 'CAC', 'CAG', 'CAT', 'AAA', 'CAA', 'GAA', 'TAA'})


if __name__ == "__main__":
    unittest.main()
