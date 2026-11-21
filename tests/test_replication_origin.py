"""Tests using the sample datasets from the Biology Meets Programming course.

Run from the project folder with:   python3 -m unittest
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from replication_origin import *  # noqa: E402,F403


class TestCounting(unittest.TestCase):
    def test_pattern_count(self):
        self.assertEqual(PatternCount("GCGCG", "GCG"), 2)
        self.assertEqual(PatternCount(ori, "TGATCA"), 8)

    def test_frequency_map(self):
        self.assertEqual(FrequencyMap("CGATATATCCATAG", 3)["ATA"], 3)

    def test_frequent_words(self):
        self.assertEqual(sorted(FrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)),
                         ["CATG", "GCAT"])
        self.assertEqual(FrequentWords("AC", 5), [])

    def test_frequent_9mers_in_ori(self):
        self.assertEqual(sorted(FrequentWords(ori, 9)),
                         ["ATGATCAAG", "CTCTTGATC", "CTTGATCAT", "TCTTGATCA"])


class TestStrings(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(ReverseComplement("AAAACCCGGT"), "ACCGGGTTTT")
        self.assertEqual(ReverseComplement("atgatcaag"), "CTTGATCAT")

    def test_pattern_matching(self):
        self.assertEqual(PatternMatching("ATAT", "GATATATGCATATACTT"), [1, 3, 9])


class TestSkew(unittest.TestCase):
    def test_skew_array(self):
        self.assertEqual(SkewArray("CATGGGCATCGGCCATACGCC"),
                         [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2])

    def test_minimum_skew(self):
        genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
        self.assertEqual(MinimumSkew(genome), [11, 24])

    def test_symbol_array(self):
        self.assertEqual(SymbolArray("AAAAGGGG", "A"),
                         {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 1, 6: 2, 7: 3})


class TestMismatches(unittest.TestCase):
    def test_hamming_distance(self):
        self.assertEqual(HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC"), 3)
        with self.assertRaises(ValueError):
            HammingDistance("ACG", "AC")

    def test_approximate_pattern_matching(self):
        text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
        self.assertEqual(ApproximatePatternMatching(text, "ATTCTGGA", 3), [6, 7, 26, 27])

    def test_approximate_pattern_count(self):
        self.assertEqual(ApproximatePatternCount("GAGG", "TTTAGAGCCTTCAGAGG", 2), 4)


class TestClumps(unittest.TestCase):
    def test_find_clumps(self):
        text = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
        self.assertEqual(FindClumps(text, 5, 50, 4), ["CGACA", "GAAGA"])

    def test_find_clumps_edge_cases(self):
        self.assertEqual(FindClumps("ACGT", 2, 10, 1), [])


if __name__ == "__main__":
    unittest.main()
