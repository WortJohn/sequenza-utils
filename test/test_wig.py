import unittest
from sequenza import wig, misc
try:
    from cStringIO import StringIO
except:
    try:
        from StringIO import StringIO
    except:
        from io import StringIO


class TestWig(unittest.TestCase):

    def test_wiggle_out(self):
        with misc.xopen('test/data/gc/summary_gc50.wig.gz', 'rt') as gc:
            wiggle = wig.Wiggle(gc)
            for line in wiggle:
                last = line
        self.assertEqual(last, (('17', 49951, 50001), ('58',)))

    def test_wiggle_in(self):
        with misc.xopen('-', 'wt') as gc:
            wiggle = wig.Wiggle(gc)
            with self.assertRaises(wig.WiggleError):
                wiggle.write((1, 10, 5, 50))
        out = StringIO()
        wiggle = wig.Wiggle(out)
        wiggle.write((1, 1, 11, 1))
        wiggle.write((1, 11, 21, 2))
        wiggle.write((2, 1, 11, 3))
        wiggle.write((2, 11, 21, 4))
        lines = StringIO(out.getvalue())
        self.assertEqual(next(lines).strip(), 'variableStep chrom=1 span=10')
        self.assertEqual(next(lines).strip().split('\t'), ['1', '1'])
        self.assertEqual(next(lines).strip().split('\t'), ['11', '2'])
        self.assertEqual(next(lines).strip(), 'variableStep chrom=2 span=10')
        self.assertEqual(next(lines).strip().split('\t'), ['1', '3'])
        self.assertEqual(next(lines).strip().split('\t'), ['11', '4'])
        out.close()

if __name__ == '__main__':
    unittest.main()
