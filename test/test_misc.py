import unittest
import argparse
from sequenza import misc
from sequenza import programs


class TestMiscs(unittest.TestCase):

    def test_xopen(self):
        file = "test/data/gc/summary_gc50.wig.gz"
        with misc.xopen(file, "rb") as gc:
            head = next(gc).strip()
        self.assertEqual(head, b'variableStep chrom=7 span=50')

    def test_countN(self):
        self.assertEqual(misc.countN('CCCCC', 'C'), 100)
        self.assertEqual(misc.countN('CCCCC', 'G'), 0)
        seq = 'AAACCCGGG'
        self.assertEqual(
            int(misc.countN(seq, 'C') + misc.countN(seq, 'G')), 66)

    def setUp(self):
        self.parser = argparse.ArgumentParser(prog='sqz', description='Test')
        self.subparsers = self.parser.add_subparsers(dest='modules')

    def test_package_modules(self):
        mods = misc.package_modules(programs)
        self.assertEqual(len(mods), 6)
        self.assertEqual(type(mods), set)
        self.assertTrue('sequenza.programs.bam2seqz' in mods)

    def test_get_modules(self):
        mods = misc.get_modules(programs, self.subparsers, {})
        self.assertEqual(hasattr(mods['bam2seqz'], '__call__'), True)

if __name__ == '__main__':
    unittest.main()
