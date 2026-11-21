"""Tests using the sample datasets from the Biology Meets Programming course,
plus edge cases, randomized checks against brute force, and menu robustness.

Run from the project folder with:   python3 -m unittest
"""
import itertools
import os
import random
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
from replication_origin import (  # noqa: E402
    ApproximatePatternCount, ApproximatePatternMatching, FindClumps, FrequencyMap,
    FrequentWords, FrequentWordsWithMismatches, FrequentWordsWithMismatchesAndReverseComplements,
    FrequentWordsWithMismatchesScores, HammingDistance, MinimumSkew, Neighbors, NeighborhoodSize,
    PatternCount, PatternMatching, ReverseComplement, SkewArray, SymbolArray, dnaa_candidates,
    ecoli_window, ori,
)


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


class TestFrequentWordsWithMismatches(unittest.TestCase):
    def test_neighbors(self):
        self.assertEqual(sorted(Neighbors("ACG", 1)),
                         ["AAG", "ACA", "ACC", "ACG", "ACT", "AGG", "ATG", "CCG", "GCG", "TCG"])
        self.assertEqual(Neighbors("ACG", 0), {"ACG"})

    def test_frequent_words_with_mismatches(self):
        text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
        self.assertEqual(FrequentWordsWithMismatches(text, 4, 1), ["ATGC", "ATGT", "GATG"])
        self.assertEqual(FrequentWordsWithMismatchesAndReverseComplements(text, 4, 1),
                         ["ACAT", "ATGT"])

    def test_scores_match_brute_force(self):
        scores = FrequentWordsWithMismatchesScores(ecoli_window, 9, 1)
        for kmer in ["TTATCCACA", "GGATCCTGG", "AAAAAAAAA", "ATGATCAAG"]:
            expected = (ApproximatePatternCount(kmer, ecoli_window, 1)
                        + ApproximatePatternCount(ReverseComplement(kmer), ecoli_window, 1))
            self.assertEqual(scores.get(kmer, 0), expected)

    def test_ecoli_dnaa_box(self):
        # Textbook: no 9-mer appears 3+ times exactly (with reverse complement)...
        exact = [k for k in FrequencyMap(ecoli_window, 9)
                 if PatternCount(ecoli_window, k)
                 + PatternCount(ecoli_window, ReverseComplement(k)) >= 3]
        self.assertEqual(exact, [])
        # ...but with 1 mismatch the real DnaA box is a most frequent 9-mer,
        # alongside the five others the textbook lists.
        best, tied, observed = dnaa_candidates(ecoli_window, 9, 1)
        self.assertEqual(best, 4)
        for kmer in ["TTATCCACA", "TGTGGATAA", "GGATCCTGG", "GATCCCAGC",
                     "GTTATCCAC", "AGCTGGGAT", "CTGGGATCA"]:
            self.assertIn(kmer, observed)

    def test_vibrio_approximate_boxes(self):
        hits = (ApproximatePatternCount("ATGATCAAG", ori, 1)
                + ApproximatePatternCount("CTTGATCAT", ori, 1))
        self.assertEqual(hits, 9)   # 8 sites; one overlapping site matches both


class TestClumps(unittest.TestCase):
    def test_find_clumps(self):
        text = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
        self.assertEqual(FindClumps(text, 5, 50, 4), ["CGACA", "GAAGA"])

    def test_find_clumps_edge_cases(self):
        self.assertEqual(FindClumps("ACGT", 2, 10, 1), [])


class TestLowercase(unittest.TestCase):
    """Lowercase DNA (common in FASTA files) used to give silently wrong answers."""

    def test_counting_and_matching(self):
        self.assertEqual(PatternCount("acgacg", "ACG"), 2)
        self.assertEqual(PatternMatching("acg", "ACGACG"), [0, 3])
        self.assertEqual(FrequentWords("aaat", 2), ["AA"])
        self.assertEqual(FrequencyMap("acac", 2), {"AC": 2, "CA": 1})

    def test_mismatch_functions(self):
        self.assertEqual(HammingDistance("acg", "ACG"), 0)
        self.assertEqual(ApproximatePatternCount("ACG", "acgacg", 0), 2)
        self.assertEqual(ApproximatePatternMatching("acgacg", "ACG", 0), [0, 3])
        self.assertEqual(Neighbors("ac", 0), {"AC"})
        self.assertIn("ACG", Neighbors("acg", 1))

    def test_skew_and_clumps(self):
        self.assertEqual(SymbolArray("ccaa", "c")[0], 2)
        self.assertEqual(MinimumSkew("cc"), [2])
        self.assertEqual(FindClumps("acacac", 2, 6, 3), ["AC"])
        self.assertEqual(ReverseComplement("acg"), "CGT")


class TestInvalidArguments(unittest.TestCase):
    """Nonsense arguments raise ValueError instead of returning an answer."""

    def test_empty_pattern(self):
        for call in (lambda: PatternCount("ACGT", ""), lambda: PatternMatching("", "ACGT"),
                     lambda: ApproximatePatternCount("", "ACGT", 1),
                     lambda: ApproximatePatternMatching("ACGT", "", 1)):
            with self.assertRaises(ValueError):
                call()

    def test_bad_k_d_t(self):
        for call in (lambda: FrequentWords("ACGT", 0), lambda: FrequencyMap("ACGT", -1),
                     lambda: FindClumps("ACGT", 0, 4, 1), lambda: FindClumps("ACGT", 2, 4, 0),
                     lambda: Neighbors("ACG", -1), lambda: ApproximatePatternCount("A", "AAA", -1),
                     lambda: FrequentWordsWithMismatches("ACGT", 0, 1)):
            with self.assertRaises(ValueError):
                call()


class TestEdgeCases(unittest.TestCase):
    def test_empty_inputs(self):
        self.assertEqual(SkewArray(""), [0])
        self.assertEqual(MinimumSkew(""), [0])
        self.assertEqual(SymbolArray("", "C"), {})
        self.assertEqual(Neighbors("", 1), {""})
        self.assertEqual(FrequentWordsWithMismatches("", 3, 1), [])

    def test_long_pattern_neighbors(self):
        # used to crash with RecursionError past ~1000 bases
        pattern = "ACGT" * 500
        self.assertEqual(len(Neighbors(pattern, 1)), NeighborhoodSize(2000, 1))

    def test_neighborhood_size(self):
        self.assertEqual(NeighborhoodSize(9, 1), 28)
        self.assertEqual(NeighborhoodSize(9, 3), 2620)
        self.assertEqual(NeighborhoodSize(2, 5), 16)

    def test_other_characters(self):
        self.assertEqual(ReverseComplement("ANT"), "ANT")
        self.assertEqual(SkewArray("GNC"), [0, 1, 1, 0])


class TestAgainstBruteForce(unittest.TestCase):
    """Random small inputs checked against slow, obviously-correct versions."""

    def setUp(self):
        self.rng = random.Random(2026)

    def dna(self, n, alphabet="ACGT"):
        return "".join(self.rng.choice(alphabet) for _ in range(n))

    @staticmethod
    def hamming(p, q):
        return sum(a != b for a, b in zip(p, q))

    def test_symbol_array_circular(self):
        for _ in range(200):
            g, s = self.dna(self.rng.randint(1, 30)), self.rng.choice("ACGT")
            n = len(g)
            expected = {i: (g + g)[i:i + n // 2].count(s) for i in range(n)}
            self.assertEqual(SymbolArray(g, s), expected)

    def test_neighbors(self):
        for _ in range(100):
            p, d = self.dna(self.rng.randint(1, 5)), self.rng.randint(0, 3)
            expected = {"".join(x) for x in itertools.product("ACGT", repeat=len(p))
                        if self.hamming(p, "".join(x)) <= d}
            self.assertEqual(Neighbors(p, d), expected)

    def test_frequent_words_with_mismatches(self):
        for _ in range(30):
            text, k, d = self.dna(self.rng.randint(4, 20)), self.rng.randint(2, 4), self.rng.randint(0, 2)
            for use_rc in (False, True):
                scores = {}
                for x in itertools.product("ACGT", repeat=k):
                    p = "".join(x)
                    s = ApproximatePatternCount(p, text, d)
                    if use_rc:
                        s += ApproximatePatternCount(ReverseComplement(p), text, d)
                    scores[p] = s
                best = max(scores.values())
                expected = sorted(p for p in scores if scores[p] == best)
                f = FrequentWordsWithMismatchesAndReverseComplements if use_rc else FrequentWordsWithMismatches
                self.assertEqual(f(text, k, d), expected)

    def test_find_clumps(self):
        for _ in range(100):
            text = self.dna(self.rng.randint(3, 40), "ACG")
            k, L, t = self.rng.randint(1, 3), self.rng.randint(3, len(text)), self.rng.randint(1, 4)
            expected = set()
            for s in range(len(text) - L + 1):
                w = text[s:s + L]
                for i in range(L - k + 1):
                    if PatternCount(w, w[i:i + k]) >= t:
                        expected.add(w[i:i + k])
            self.assertEqual(FindClumps(text, k, L, t), sorted(expected) if k <= L else [])


class TestMenu(unittest.TestCase):
    """Drive the interactive menu with typed input; it should never crash or hang."""

    def run_menu(self, keys, files=None):
        with tempfile.TemporaryDirectory() as work:
            for name, content in (files or {}).items():
                with open(os.path.join(work, name), "wb") as f:
                    f.write(content)
            script = os.path.join(ROOT, "replication_origin.py")
            result = subprocess.run([sys.executable, script], input="\n".join(keys) + "\n",
                                    capture_output=True, text=True, timeout=60, cwd=work)
        self.assertNotIn("Traceback", result.stdout + result.stderr)
        return result.stdout

    def test_every_option(self):
        out = self.run_menu(["1", "1", "ATG", "2", "9", "3", "acgt", "4", "ATGATCAAG", "5", "C", "6", "7",
                             "9", "ACG", "ACT", "10", "ATGATCAAG", "1", "11", "ATGATCAAG", "1",
                             "12", "9", "1", "s", "13", "9", "500", "3", "14", "9", "1", "y", "w",
                             "s", "1", "q"])
        self.assertTrue(out.rstrip().endswith("Goodbye!"))

    def test_bad_typing(self):
        out = self.run_menu(["1", "99", "", "2", "abc", "-3", "99999", "9", "9", "ACGT", "AC", "q"])
        self.assertIn("Not a valid option", out)
        self.assertIn("same length", out)

    def test_bad_files(self):
        out = self.run_menu(["6", "pic.png", "s", "6", "empty.txt", "s", "6", "missing.txt", "q"],
                            files={"pic.png": bytes(range(256)) * 10, "empty.txt": b""})
        self.assertIn("doesn't look like a text file", out)
        self.assertIn("has no sequence", out)
        self.assertIn("Could not open", out)

    def test_huge_mismatch_search_asks_first(self):
        out = self.run_menu(["1", "14", "9", "9", "y", "n", "q"])
        self.assertIn("Heads up", out)
        self.assertIn("Cancelled", out)

    def test_ctrl_d_quits_cleanly(self):
        self.assertIn("Goodbye", self.run_menu(["1", "2"]))


if __name__ == "__main__":
    unittest.main()
