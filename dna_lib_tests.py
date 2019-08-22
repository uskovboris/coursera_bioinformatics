import unittest
from dna_lib import *


class CheckNucleotideTest(unittest.TestCase):

    def test_check_nucleotide_adenine_ok(self):
        """ Check A """
        self.assertTrue(check_nucleotide("A"), "Adenine is nucleotide")

    def test_check_nucleotide_timin_ok(self):
        """ Check T """
        self.assertTrue(check_nucleotide("T"), "Timin is nucleotide")

    def test_check_nucleotide_cytosine_ok(self):
        """ Check C """
        self.assertTrue(check_nucleotide("C"), "Cytosine is nucleotide")

    def test_check_nucleotide_guanine_ok(self):
        """ Check G """
        self.assertTrue(check_nucleotide("G"), "Guanine is nucleotide")

    def test_check_nucleotide_guanine_ok(self):
        """ Check not-existed nucleotide """
        self.assertFalse(check_nucleotide("Y"), "Y is not nucleotide")


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

    def test_hamming_dist_strings(self):
        self.assertEqual(hamming_dist("GGGCCGTTGGT", "GGACCGTTGAC"), 3)


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


class IterativeNeighboursTest(unittest.TestCase):

    def test_neighbours_d_equals_0(self):
        """d=0 test"""
        self.assertSetEqual(iterative_neighbors('ACG', 0), {"ACG"})

    def test_neighbours_len_equals_1(self):
        """1 nucleotide len DNA test"""
        self.assertSetEqual(iterative_neighbors('A', 1), {"A", "T", "C", "G"})

    def test_neighbours(self):
        """ neighbours test """
        self.assertSetEqual(iterative_neighbors('ACG', 1), {'CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 'ACA', 'ACC', 'ACT', 'ACG'})
        self.assertSetEqual(iterative_neighbors('CAA', 1), {'CAA', 'CCA', 'CGA', 'CTA', 'CAC', 'CAG', 'CAT', 'AAA', 'CAA', 'GAA', 'TAA'})


if __name__ == "__main__":
    unittest.main()
