import unittest
from dna_lib import *


class CheckNucleotideTest(unittest.TestCase):

    def test_check_nucleotide_adenine_ok(self):
        """ Check A """
        self.assertTrue(check_nucleotide("A"), "Adenine is nucleotide")

    def test_check_nucleotide_thymine_ok(self):
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


class SkewTest(unittest.TestCase):

    def test_skew_i_empty(self):
        with self.assertRaises(ValueError):
            skew("", 0)

    def test_skew_i_less_then_0(self):
        with self.assertRaises(ValueError):
            skew("ATCG", -1)

    def test_skew_i_equal_0(self):
        self.assertEqual(skew("CATGGGCATCGGCCATACGCC", 0), 0)

    def test_skew_cytosine(self):
        self.assertEqual(skew("C", 1), -1)

    def test_skew_guanine(self):
        self.assertEqual(skew("G", 1), 1)

    def test_skew_adenine(self):
        self.assertEqual(skew("GA", 2), 1)

    def test_skew_thymine(self):
        self.assertEqual(skew("GT", 2), 1)

    def test_check_nucleotide_adenine_ok(self):
        """ Complex test"""
        pattern = "CATGGGCATCGGCCATACGCC"
        self.assertEqual(skew(pattern, 0), 0)
        self.assertEqual(skew(pattern, 1), -1)
        self.assertEqual(skew(pattern, 2), -1)
        self.assertEqual(skew(pattern, 3), -1)
        self.assertEqual(skew(pattern, 4), 0)
        self.assertEqual(skew(pattern, 5), 1)
        self.assertEqual(skew(pattern, 6), 2)
        self.assertEqual(skew(pattern, 7), 1)
        self.assertEqual(skew(pattern, 8), 1)
        self.assertEqual(skew(pattern, 9), 1)


class FindSkewMinimumsTest(unittest.TestCase):

    def test_genome_is_empty(self):
        with self.assertRaises(ValueError):
            find_skew_minimums("")

    def test_find_skew_minimum_dataset1(self):
        """
        The sample dataset is not actually run on your code
        """
        indexes, _ = find_skew_minimums("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")
        self.assertSequenceEqual([11, 24], indexes)

    def test_find_skew_minimum_indexing_is_off(self):
        """
        This dataset checks if your code’s indexing is off. Specifically, it verifies that your code
        is not returning an index 1 too high (i.e. 4) or 1 too low (i.e. 2).
        """
        indexes, _ = find_skew_minimums("ACCG")
        self.assertSequenceEqual([3], indexes)

    def test_find_skew_minimum_check_last_symbol_in_genome(self):
        """
        This dataset checks to see if your code is missing the last symbol of Genome.
        """
        indexes, _ = find_skew_minimums("ACCC")
        self.assertSequenceEqual([4], indexes)

    def test_find_skew_minimum_check_there_is_minimum_not_maximum(self):
        """
        This dataset makes sure you’re not accidentally finding the maximum skew instead of the
        minimum skew
        """
        indexes, _ = find_skew_minimums("CCGGGT")
        self.assertSequenceEqual([2], indexes)


if __name__ == "__main__":
    unittest.main()
