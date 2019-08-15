import unittest
from clump_finding import *


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


if __name__ == "__main__":
    unittest.main()
