import unittest
from sequenza.misc import xopen
from sequenza.seqz import do_seqz, seqz_header


class TestSeqz(unittest.TestCase):

    # def setUp(self):
    #     with xopen('test/data/mpileup/testnorm.pileup.gz', 'rt') as test:
    #         for line in test:
    #             if line.startswith('12\t299922\t'):
    #                 self.normal = line
    #     with xopen('test/data/mpileup/testtum.pileup.gz', 'rt') as test:
    #         for line in test:
    #             if line.startswith('12\t299922\t'):
    #                 self.tumor = line
    #     self.normal = '\t'.join(self.normal.rstrip().split('\t')[2:])
    #     self.tumor = '\t'.join(self.tumor.rstrip().split('\t')[2:])

    def test_do_seqz(self):
        self.normal = ('T\t29\t,C.C,c,C,c,,,c,cCccC,c,,c,c,,\tBB/'
                       '<FFFBFFFFFFBFFFF/7/7FBFFFF')
        self.tumor = ('T\t46\tc$ccc,cCcc.cc,cGcC.Ccc,c,C.CC,.CcC.ccc,Cc,ccccg'
                      '\t/FFFFFgBF/FF/F/F<//FF/FFk</BF/<F/BF/FFB/FBFF</')
        data = (self.normal, self.tumor, '50')
        seqz_list = do_seqz(data)
        self.assertEqual(seqz_list[9], 'CT')
        self.normal = ('T\t29\t,...,,,.,,,,,,,,.,,.,,,,,,,,,\tBB/'
                       '<FFFBFFFFFFBFFFF/7/7FBFFFF')
        self.tumor = ('T\t46\tc$ccc,cCcc.cc,cGcC.Ccc,c,C.CC,.CcC.ccc,Cc,ccccg'
                      '\t/FFFFFgBF/FF/F/F<//FF/FFk</BF/<F/BF/FFB/FBFF</')
        data = (self.normal, self.tumor, '50')
        seqz_list = do_seqz(data)
        self.assertEqual(seqz_list[9], 'T')
        self.assertEqual(seqz_list[10], 'C0.788')
        self.assertEqual(seqz_list[11], 'C0.231')
        data = (self.normal, self.tumor, '50', self.tumor)
        seqz_list = do_seqz(data)
        self.assertEqual(seqz_list[1], seqz_list[2])
        self.assertEqual(seqz_list[9], 'T')
        self.assertEqual(seqz_list[10], 'C0.788')
        self.assertEqual(seqz_list[11], 'C0.231')

    def test_header(self):
        header = seqz_header()
        self.assertEqual(len(header), 14)


if __name__ == '__main__':
    unittest.main()
