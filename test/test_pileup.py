import unittest
from sequenza.misc import xopen
from sequenza import c_pileup, pileup


class TestPileup(unittest.TestCase):

    def setUp(self):
        with xopen('test/data/mpileup/testnorm.pileup.gz', 'rt') as test:
            for line in test:
                if line.startswith('12\t240965\t'):
                    self.test_het = line
                if line.startswith('17\t41629\t'):
                    self.test_mut = line
        self.test_het = self.test_het.split('\t')
        self.test_mut = self.test_mut.split('\t')

    def test_pileup_parser(self):
        test = pileup.acgt(pileup=self.test_het[4], quality=self.test_het[5],
                           depth=int(self.test_het[3]),
                           reference=self.test_het[2])
        self.assertEqual(test['Z'], [0, 31, 43, 0])
        self.assertEqual(test['A'], 0)
        self.assertEqual(test['C'], 36)
        self.assertEqual(test['G'], 44)
        self.assertEqual(test['T'], 0)

        test = pileup.acgt(pileup=self.test_mut[4], quality=self.test_mut[5],
                           depth=int(self.test_mut[3]),
                           reference=self.test_mut[2])
        self.assertEqual(test['Z'], [0, 3, 0, 1])
        self.assertEqual(test['A'], 0)
        self.assertEqual(test['C'], 16)
        self.assertEqual(test['G'], 0)
        self.assertEqual(test['T'], 13)

        test = c_pileup.acgt(pileup=self.test_het[4], quality=self.test_het[5],
                             depth=int(self.test_het[3]),
                             reference=self.test_het[2])
        self.assertEqual(test['Z'], [0, 31, 43, 0])
        self.assertEqual(test['A'], 0)
        self.assertEqual(test['C'], 36)
        self.assertEqual(test['G'], 44)
        self.assertEqual(test['T'], 0)

        test = c_pileup.acgt(pileup=self.test_mut[4], quality=self.test_mut[5],
                             depth=int(self.test_mut[3]),
                             reference=self.test_mut[2])
        self.assertEqual(test['Z'], [0, 3, 0, 1])
        self.assertEqual(test['A'], 0)
        self.assertEqual(test['C'], 16)
        self.assertEqual(test['G'], 0)
        self.assertEqual(test['T'], 13)


if __name__ == '__main__':
    unittest.main()
