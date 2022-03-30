import random
from unittest import TestCase

from archive2dna import Dna


class DnaClass(TestCase):
    def test_init_empty(self):
        self.assertEqual(Dna().to_string(), "")
        self.assertEqual(Dna().to_bytes(), bytes())

    def test_bytes_to_string(self):
        def _test_bytes_to_string(b, d):
            self.assertEqual(Dna(b).to_string(), d)
            self.assertEqual(Dna.from_bytes(b).to_string(), d)

        _test_bytes_to_string(bytes(), "")
        _test_bytes_to_string(bytes([0b11100100]), "TCGA")
        _test_bytes_to_string(bytes([0b10000111, 0b01001110]), "CAGTGATC")

    def test_string_to_bytes(self):
        def _test_string_to_bytes(d, b):
            self.assertEqual(Dna(d).to_bytes(), b)
            self.assertEqual(Dna.from_string(d).to_bytes(), b)

        _test_string_to_bytes("", bytes())
        _test_string_to_bytes("TGAC", bytes([0b11010010]))
        _test_string_to_bytes("GGTCATCG", bytes([0b01011110, 0b00111001]))

    def test_from_string_validation(self):
        with self.assertRaisesRegex(
            Exception, "dna expected to be a multiple of 4 bases"
        ):
            Dna("C")
        with self.assertRaisesRegex(
            Exception, "dna expected to be a multiple of 4 bases"
        ):
            Dna.from_string("C")
        with self.assertRaisesRegex(
            Exception, "dna expected to be a multiple of 4 bases"
        ):
            Dna("CATGDC")
        with self.assertRaisesRegex(
            Exception, "dna expected to be a multiple of 4 bases"
        ):
            Dna.from_string("CATGDC")
        with self.assertRaisesRegex(Exception, "base expected to be A, T, G or C"):
            Dna("CATV")
        with self.assertRaisesRegex(Exception, "base expected to be A, T, G or C"):
            Dna.from_string("CATV")

    def test_bytes_to_string_and_back(self):
        b = bytes(random.getrandbits(8) for _ in range(256))
        self.assertTrue(Dna(Dna(b).to_string()).to_bytes() == b)
        d = "".join([random.choice("ATCG") for _ in range(256)])
        self.assertTrue(Dna(d).to_string() == d)

    def test_equality(self):
        self.assertEqual(Dna(""), Dna(""))
        self.assertEqual(Dna("ACTG"), Dna("ACTG"))
        self.assertNotEqual(Dna("ACTG"), Dna("ACTC"))
        self.assertEqual(Dna("ACTGTTCA"), Dna("ACTGTTCA"))
        self.assertNotEqual(Dna("ACTGTTCA"), Dna("ACTGTGCA"))

    def test_len(self):
        self.assertEqual(len(Dna("")), 0)
        self.assertEqual(len(Dna("ATCG")), 4)
        self.assertEqual(len(Dna("ATCGTCGG")), 8)
        self.assertEqual(len(Dna("ATCGTCGGCACC")), 12)

    def test_get(self):
        d = Dna("ATCGCTTA")
        self.assertEqual(d[0], "A")
        self.assertEqual(d[1], "T")
        self.assertEqual(d[2], "C")
        self.assertEqual(d[3], "G")
        self.assertEqual(d[4], "C")
        self.assertEqual(d[5], "T")
        self.assertEqual(d[6], "T")
        self.assertEqual(d[7], "A")
        with self.assertRaises(IndexError):
            d[8]

    def test_concat(self):
        def _test_concat(d1, d2, result):
            self.assertEqual(d1.concat(d2), result)
            self.assertEqual(d1 + d2, result)
            d1 += d2
            self.assertEqual(d1, result)

        _test_concat(Dna(""), Dna(""), Dna(""))
        _test_concat(Dna(""), Dna("ATGG"), Dna("ATGG"))
        _test_concat(Dna("CTAG"), Dna(""), Dna("CTAG"))
        _test_concat(Dna("CGAA"), Dna("GTCA"), Dna("CGAAGTCA"))

    def test_mul(self):
        def _test_mul(d, i, result):
            self.assertEqual(d * i, result)
            self.assertEqual(i * d, result)
            d *= i
            self.assertEqual(d, result)

        _test_mul(Dna(""), 0, Dna(""))
        _test_mul(Dna(""), 1, Dna(""))
        _test_mul(Dna(""), 2, Dna(""))
        _test_mul(Dna("CTAG"), 0, Dna(""))
        _test_mul(Dna("CTAG"), 1, Dna("CTAG"))
        _test_mul(Dna("CTAG"), 2, Dna("CTAGCTAG"))
        _test_mul(Dna("CGAATCTG"), 0, Dna(""))
        _test_mul(Dna("CGAATCTG"), 1, Dna("CGAATCTG"))
        _test_mul(Dna("CGAATCTG"), 2, Dna("CGAATCTGCGAATCTG"))

    def test_complement(self):
        def _test_complement(d, result):
            self.assertEqual(d.complement(), result)
            self.assertEqual(~d, result)

        _test_complement(Dna(""), Dna(""))
        _test_complement(Dna("CTGA"), Dna("GACT"))
        _test_complement(Dna("CGAATCTG"), Dna("GCTTAGAC"))

    def test_addPrimer(self):
        self.assertEqual(Dna("").addPrimer(Dna("")), Dna(""))
        self.assertEqual(Dna("").addPrimer(Dna("AGGT")), Dna("AGGTTCCA"))
        self.assertEqual(
            Dna("TCGTTGCA").addPrimer(Dna("CTGT")), Dna("CTGTTCGTTGCAGACA")
        )

    def test_removePrimer(self):
        self.assertEqual(Dna("").removePrimer(Dna("")), Dna(""))
        self.assertEqual(Dna("AGGTTCCA").removePrimer(Dna("AGGT")), Dna(""))
        self.assertEqual(
            Dna("CTGTTCGTTGCAGACA").removePrimer(Dna("CTGT")), Dna("TCGTTGCA")
        )
        with self.assertRaisesRegex(Exception, "start primer not detected"):
            Dna("CTATTCGTTGCAGACA").removePrimer(Dna("CTGT"))
        with self.assertRaisesRegex(Exception, "end primer not detected"):
            Dna("CTGTTCGTTGCAGACT").removePrimer(Dna("CTGT"))
