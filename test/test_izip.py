import unittest
from sequenza import wig, misc, izip


class TestIzip(unittest.TestCase):

    def test_coordinate_zip(self):
        with misc.xopen('test/data/mpileup/testnorm.pileup.gz', 'rt') \
                as normal, \
                misc.xopen('test/data/mpileup/testtum.pileup.gz', 'rt') \
                as tumor, \
                misc.xopen('test/data/gc/summary_gc50.wig.gz', 'rt') \
                as gc:
            nor_tum = izip.zip_coordinates(misc.split_coordinates(
                normal), misc.split_coordinates(tumor))
            nor_tum_gc = izip.zip_coordinates(nor_tum, wig.Wiggle(gc))
            res = nor_tum_gc.next()
            self.assertEqual(res[1][2], '60')
            self.assertEqual(res[1][0].split('\t')[0],
                             res[1][1].split('\t')[0])
            self.assertEqual(res[0][1], 223)
            res = next(nor_tum_gc)
            self.assertEqual(res[1][2], '60')
            self.assertEqual(res[1][0].split('\t')[0],
                             res[1][1].split('\t')[0])
            self.assertEqual(res[0][1], 224)
            nor_tum_gc.close()


if __name__ == '__main__':
    unittest.main()
