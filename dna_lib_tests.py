import unittest
from dna_lib import *


class NucleotideToNumberTest(unittest.TestCase):

    def test_dna_adenine_to_number(self):
        self.assertEqual(nucleotide_to_number('A'), 0)

    def test_dna_thymine_to_number(self):
        self.assertEqual(nucleotide_to_number('T'), 3)

    def test_dna_cytosine_to_number(self):
        self.assertEqual(nucleotide_to_number('C'), 1)

    def test_dna_guanine_to_number(self):
        self.assertEqual(nucleotide_to_number('G'), 2)

    def test_rna_adenine_to_number(self):
        self.assertEqual(nucleotide_to_number('A', RNA=True), 0)

    def test_rna_uracile_to_number(self):
        self.assertEqual(nucleotide_to_number('U', RNA=True), 3)

    def test_rna_cytosine_to_number(self):
        self.assertEqual(nucleotide_to_number('C', RNA=True), 1)

    def test_rna_guanine_to_number(self):
        self.assertEqual(nucleotide_to_number('G', RNA=True), 2)


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
        self.assertSequenceEqual(neighbors('ACG', 0), ["ACG"])

    def test_neighbours_len_equals_1(self):
        """1 nucleotide len DNA test"""
        self.assertSequenceEqual(sorted(neighbors('A', 1)), sorted(["A", "T", "C", "G"]))

    def test_neighbours(self):
        """ neighbours test """
        self.assertSequenceEqual(sorted(neighbors('ACG', 1)),
                                 sorted(['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 'ACA', 'ACC', 'ACT', 'ACG']))
        self.assertSequenceEqual(sorted(neighbors('CAA', 1)),
                                 sorted(['CAA', 'CCA', 'CGA', 'CTA', 'CAC', 'CAG', 'CAT', 'AAA', 'GAA', 'TAA']))


class IterativeNeighboursTest(unittest.TestCase):

    def test_neighbours_d_equals_0(self):
        """d=0 test"""
        self.assertSequenceEqual(sorted(iterative_neighbors('ACG', 0)), sorted(["ACG"]))

    def test_neighbours_len_equals_1(self):
        """1 nucleotide len DNA test"""
        self.assertSequenceEqual(sorted(iterative_neighbors('A', 1)), sorted(["A", "T", "C", "G"]))

    def test_neighbours(self):
        """ neighbours test """
        self.assertSequenceEqual(sorted(iterative_neighbors('ACG', 1)),
                                 sorted(['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 'ACA', 'ACC', 'ACT', 'ACG']))
        self.assertSequenceEqual(sorted(iterative_neighbors('CAA', 1)),
                                 sorted(['CAA', 'CCA', 'CGA', 'CTA', 'CAC', 'CAG', 'CAT', 'AAA', 'GAA', 'TAA']))


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


class ApproximatePatternsCountTest(unittest.TestCase):

    def test_approximate_patterns_count(self):
        """
        The sample dataset is not actually run on your code.
        """
        self.assertEqual(approximate_patterns_count("TTTAGAGCCTTCAGAGG", "GAGG", 2), 4)

    def test_approximate_patterns_overlapping(self):
        """
        Checks that function is correctly handling overlapping occurrences (i.e. returning 2 instead of 1).
        """
        self.assertEqual(approximate_patterns_count("AAA", "AA", 0), 2)

    def test_approximate_patterns_handling_patterns_with_less_then_d_mismatches(self):
        """
        This dataset checks if your code is allowing occurrences with less than d mismatches
        (which it should). It is a common mistake to only allow occurrences with exactly d mismatches,
        whereas we want all occurrences with less than or equal to d mismatches.
        """
        self.assertEqual(approximate_patterns_count("ATA", "ATA", 1), 1)

    def test_approximate_patterns_(self):
        self.assertEqual(approximate_patterns_count("ACGTTGCATGTCGCATGATGCATGAGAGCT", "ATA", 1), 1)



class PatternPositionsTest(unittest.TestCase):

    def test_patterns_positions_check_several_occurrences(self):
        self.assertSequenceEqual([1, 3, 9], pattern_positions("ATAT", "GATATATGCATATACTT"))

    def test_patterns_positions_sample2(self):
        self.assertSequenceEqual([4], pattern_positions("ACAC", "TTTTACACTTTTTTGTGTAAAAA"))

    def test_patterns_positions_check_at_begin(self):
        self.assertSequenceEqual([0, 46, 51, 74], pattern_positions("AAA",
                                                                    "AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTACACAACATCCAT"))

    def test_patterns_positions_check_tail(self):
        self.assertSequenceEqual([88, 92, 98, 132], pattern_positions("TTT",
                                                                      "AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT"))

    def test_patterns_positions_check_overlapping(self):
        self.assertSequenceEqual([0, 2, 4], pattern_positions("ATA", "ATATATA"))


class ApproximatePatternPositionsTest(unittest.TestCase):

    def test_approximate_pattern_positions(self):
        """
        The sample dataset is not actually run on your code
        """
        self.assertSequenceEqual([6, 7, 26, 27], approximate_pattern_positions("ATTCTGGA",
                                                                               "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT",
                                                                               3))

    def test_approximate_patterns_handling_patterns_with_exactly_equals_d_mismatches(self):
        """
        This dataset checks if you are only counting instances where the number of mismatches is
        exactly equal to d (i.e. ignoring instances where mismatch < d).
        """
        self.assertSequenceEqual([4, 5, 6, 7, 8, 11, 12, 13, 14, 15],
                                 approximate_pattern_positions("AAA", "TTTTTTAAATTTTAAATTTTTT", 2))

    def test_approximate_patterns_handling_patterns_check_begin(self):
        """
        This dataset checks if your code has an off-by-one error at the beginning of Text (i.e. your
        code is not checking the the leftmost substring of Text).
        """
        self.assertSequenceEqual([0, 30, 66],
                                 approximate_pattern_positions("GAGCGCTGG",
                                                               "GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGAAACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCCGCTATCGCATTTTGATGCGCGCCGCGTCGATT",
                                                               2))

    def test_approximate_patterns_handling_patterns_check_end(self):
        """
        This dataset checks if your code has an off-by-one error at the end of Text (i.e. your code
        is not checking the the rightmost substring of Text).
        """
        self.assertSequenceEqual([3, 36, 74, 137],
                                 approximate_pattern_positions("AATCCTTTCA",
                                                               "CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTCTGCATATAAGTAGTGAAGGTATAGAAACCCGTTCAAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGCGCCATAATCCAAACA",
                                                               3))

    def test_approximate_patterns_handling_patterns_check_overlapping(self):
        """
        This dataset checks if your code is correctly accounting for overlapping instances of
        Pattern in Text.
        """

        self.assertSequenceEqual([0, 7, 36, 44, 48, 72, 79, 112],
                                 approximate_pattern_positions("CCGTCATCC",
                                                               "CCGTCATCCGTCATCCTCGCCACGTTGGCATGCATTCCGTCATCCCGTCAGGCATACTTCTGCATATAAGTACAAACATCCGTCATGTCAAAGGGAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGC",
                                                               3))

    def test_approximate_patterns_handling_patterns_with_less_then_d_mismatches(self):
        """
        This dataset checks if you are only counting instances of Pattern with less than d
        mismatches (as opposed to instances of Pattern with less than or equal to d mismatches).
        """
        self.assertSequenceEqual([0, 1, 2, 3],
                                 approximate_pattern_positions("TTT",
                                                               "AAAAAA",
                                                               3))

    def test_approximate_patterns_handling_patterns_with_d_equals_0(self):
        """
        This dataset checks if your code works with input where d = 0 (i.e. only perfect matches
        are allowed).
        """
        self.assertSequenceEqual([0],
                                 approximate_pattern_positions("CCA",
                                                               "CCACCT",
                                                               0))


class ComputingFrequentPatternsWithMismatchesTest(unittest.TestCase):

    def test_computing_frequent_patterns_with_mismatches_case1(self):
        self.assertSequenceEqual(sorted(["ATGC", "ATGT", "GATG"]),
                                 sorted(computing_frequent_patterns_with_mismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)))

    def test_computing_frequent_patterns_with_mismatches_check_kmers_not_Presented_in_text(self):
        """
        This dataset checks that your code includes k-mers that do not actually appear in Text.
        Notice here that, although none of the output k-mers except for AA actually appear in Text, they
        are all valid because they appear in Text with up to 1 mismatch (i.e. 0 or 1 mismatch).
        """
        self.assertSequenceEqual(sorted(["AA", "AC", "AG", "CA", "AT", "GA", "TA"]),
                                 sorted(computing_frequent_patterns_with_mismatches("AAAAAAAAAA", 2, 1)))

    def test_computing_frequent_patterns_with_mismatches_check_k_d_positions(self):
        """
        This dataset makes sure that your code is not accidentally swapping k and d.
        """
        self.assertSequenceEqual(sorted(['AAGC', 'AATT', 'ACAC', 'ACTG', 'AGAG', 'AGCA', 'AGGT', 'ATCC', 'ATTA', 'CATC', 'CGGC', 'CGTT', 'GGCC', 'GGTA', 'GTTC', 'TCTC', 'TGAC', 'TGTG']),
                                 sorted(computing_frequent_patterns_with_mismatches("AGTCAGTC", 4, 2)))

    def test_computing_frequent_patterns_with_mismatches_check_that_not_in_reverse_complement(self):
        """
        This dataset makes sure you are not finding patterns in the Reverse Complement of Text
        """
        self.assertSequenceEqual(sorted(["GGTA"]),
                                 sorted(computing_frequent_patterns_with_mismatches("AATTAATTGGTAGGTAGGTA", 4, 0)))

    def test_computing_frequent_patterns_with_mismatches_check_that_patterns_with_less_or_equals_d_mismatches(self):
        """
        This dataset first checks that k-mers with exactly d mismatches are being found. Then, it
        checks that k-mers with less than d mismatches are being allowed (i.e. you are not only allowing
        k-mers with exactly d mismatches). Next, it checks that you are not returning too few k-mers.
        Last, it checks that you are not returning too many k-mers
        """
        self.assertSequenceEqual(sorted(["GTA", "ACA", "AAA", "ATC", "ATA", "AGA", "ATT", "CTA", "TTA", "ATG"]),
                                 sorted(computing_frequent_patterns_with_mismatches("ATA", 3, 1)))

    def test_computing_frequent_patterns_with_mismatches_check_not_in_reverse_complement(self):
        """
        This dataset checks that your code is not looking for k-mers in the Reverse Complement of Text.
        """
        self.assertSequenceEqual(sorted(["AAT"]),
                                 sorted(computing_frequent_patterns_with_mismatches("AAT", 3, 0)))


if __name__ == "__main__":
    unittest.main()
