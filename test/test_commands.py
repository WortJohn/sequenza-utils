import unittest
import sys
import os
import sequenza.commands


class TestCommands(unittest.TestCase):

    def setUp(self):
        self.vcf = 'test/data/mpileup/test.vcf.gz'
        self.bamn = 'test/data/bam/testnorm.bam'
        self.bamt = 'test/data/bam/testtum.bam'
        self.mpn = 'test/data/mpileup/testnorm.pileup.gz'
        self.mpt = 'test/data/mpileup/testnorm.pileup.gz'
        self.fasta = 'test/data/fasta/subset.fa.gz'
        self.pregc50 = 'test/data/fasta/subset.gc50.gz'

        self.gc50 = 'test/data/tmp/test_gc.wig.gz'
        self.seqz_mpi = 'test/data/tmp/test_mpi.seqz.gz'
        self.seqz_bam = 'test/data/tmp/test_bam.seqz.gz'

        self.acgt = 'test/data/tmp/test_acgt.gz'

        self.seqz_vcf = 'test/data/tmp/test_vcf.seqz.gz'
        self.seqz_vcf_bin = 'test/data/tmp/test_vcf_bin200.seqz.gz'
        self.argv = sys.argv

    def test_gc(self):
        sys.argv = ['commands', 'gc_wiggle', '-f',
                    self.fasta, '-w', '50', '-o', self.gc50]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.gc50), True)

    def test_bam_seqz(self):
        sys.argv = ['commands', 'bam2seqz', '-F', self.fasta,
                    '-gc', self.pregc50, '-t', self.bamt, '-n', self.bamn,
                    '-o', self.seqz_bam]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.seqz_bam), True)
        self.assertEqual(os.path.isfile('%s.tbi' % self.seqz_bam), True)

    def test_bam_seqz1(self):
        part_seqz = self.seqz_bam.replace('.seqz.gz', 'part_12.seqz.gz')
        sys.argv = ['commands', 'bam2seqz', '-F', self.fasta,
                    '-gc', self.pregc50, '-t', self.bamt, '-n', self.bamn,
                    '-C', '12', '-o', part_seqz]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.seqz_bam), True)
        self.assertEqual(os.path.isfile('%s.tbi' % self.seqz_bam), True)
        '''
        sys.argv = ['commands', 'bam2seqz', '-F', self.fasta,
                    '-gc', self.pregc50, '-t', self.bamt, '-n', self.bamn,
                    '-C', '7', '12', '--parallel', '2', '-o', self.seqz_bam]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.seqz_bam.replace(
            '.seqz.gz', '_7.seqz.gz')), True)
        self.assertEqual(os.path.isfile(self.seqz_bam.replace(
            '.seqz.gz', '_12.seqz.gz')), True)
        '''

    def test_mpileup_seqz(self):
        sys.argv = ['commands', 'bam2seqz', '--pileup',
                    '-gc', self.pregc50, '-t', self.mpt, '-n', self.mpn,
                    '-o', self.seqz_mpi]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.seqz_mpi), True)
        self.assertEqual(os.path.isfile('%s.tbi' % self.seqz_mpi), True)
        part_seqz = self.seqz_mpi.replace('.seqz.gz', 'part_12.seqz.gz')
        sys.argv = ['commands', 'bam2seqz', '--pileup',
                    '-gc', self.pregc50, '-t', self.mpt, '-n', self.mpn,
                    '-C', '12', '-o', part_seqz]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.seqz_mpi), True)
        self.assertEqual(os.path.isfile('%s.tbi' % self.seqz_mpi), True)
        '''
        sys.argv = ['commands', 'bam2seqz', '--pileup',
                    '-gc', self.pregc50, '-t', self.mpt, '-n', self.mpn,
                    '-C', '7', '12', '--parallel', '2', '-o', self.seqz_mpi]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.seqz_mpi.replace(
            '.seqz.gz', '_7.seqz.gz')), True)
        self.assertEqual(os.path.isfile(self.seqz_mpi.replace(
            '.seqz.gz', '_12.seqz.gz')), True)
        '''

    def test_acgt(self):
        sys.argv = ['commands', 'pileup2acgt', '--mpileup', self.mpt,
                    '-o', self.acgt]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.acgt), True)

    def test_snp(self):
        sys.argv = ['commands', 'snp2seqz', '-v', self.vcf,
                    '-gc', self.pregc50, '-o', self.seqz_vcf]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.seqz_vcf), True)
        self.assertEqual(os.path.isfile('%s.tbi' % self.seqz_vcf), True)

    def test_vcf_bin(self):
        sys.argv = ['commands', 'seqz_binning', '-w', '200',
                    '-s', self.seqz_vcf, '-o', self.seqz_vcf_bin]
        sequenza.commands.main()
        self.assertEqual(os.path.isfile(self.seqz_vcf_bin), True)
        self.assertEqual(os.path.isfile('%s.tbi' % self.seqz_vcf_bin), True)

    def tearDown(self):
        sys.argv = self.argv


if __name__ == '__main__':
    unittest.main()
