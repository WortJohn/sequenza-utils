import unittest
from sequenza import fasta, misc


class TestFasta(unittest.TestCase):

    def test_fasta_read(self):
        with misc.xopen('test/data/fasta/subset.fa.gz', 'rt') as fa:
            test = fasta.Fasta(fa, 6)
            for line in test:
                if line:
                    self.assertEqual(len(line), 4)
                last = line
            self.assertEqual(last, ('17', 49999, 50001, 'CTG'))
        with misc.xopen('test/data/fasta/subset.fa.gz', 'rt') as fa:
            test = fasta.Fasta(fa, 50)
            for line in test:
                if line:
                    self.assertEqual(len(line), 4)
                last = line
            self.assertEqual(last, ('17', 50001, 50001, 'G'))
            test.close()

if __name__ == '__main__':
    unittest.main()
