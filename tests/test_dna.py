from unittest import TestCase

from archive2dna import dna


class DnaModule(TestCase):
    def test_valid_dna_1(self):
        self.assertTrue(dna.isValidDna("ATGCATGC"))

    def test_valid_dna_2(self):
        self.assertFalse(dna.isValidDna("ATGCXTGC"))

    def test_primer_complement(self):
        self.assertTrue(dna.complement_primer("ATGC") == "GCAT")

    def test_bytes_to_dna_and_back(self):
        """Runs bytes2dna and then dna2bytes on all 256 bytes and checks
        that the result is identity (i.e. the same as the input)"""
        b = b""
        for i in range(256):
            b += dna.int2bytes(i, n=1)
        d = dna.bytes2dna(b)
        self.assertTrue(dna.dna2bytes(d) == b)

    def test_add_remove_primer(self):
        sequence = "ATGC"
        primer1 = "AAAAAA"
        primer2 = "CCCCCC"
        s2 = dna.add_primers(sequence, primer1=primer1, primer2=primer2)
        s3 = dna.remove_primers(s2, primer1=primer1, primer2=primer2)
        self.assertTrue(sequence == s3)
