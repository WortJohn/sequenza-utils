###############################
#
# Prepare test files:
#  - Use bwa to build BAM files
#  - Use samtools to build pileups
#  - oters
###############################

# External program requirements:
# BWA mem
# samtools
# tabix
# bgzip

import os
from shutil import rmtree
import shlex
import subprocess
from sequenza.samtools import bam_mpileup
from sequenza.misc import xopen, countN
from sequenza.fasta import Fasta
from sequenza.wig import Wiggle

# Set path and file names
test_data = 'test/data'
tmp_dir = os.path.join(test_data, 'tmp')
fastq_dir = os.path.join(test_data, 'fastq')
fasta_dir = os.path.join(test_data, 'fasta')
bam_dir = os.path.join(test_data, 'bam')
pup_dir = os.path.join(test_data, 'mpileup')
gc_dir = os.path.join(test_data, 'gc')

ref_fasta = 'subset.fa.gz'
bam_tumor = 'testtum.bam'
bam_normal = 'testnorm.bam'
pup_tumor = 'testtum.pileup'
pup_normal = 'testnorm.pileup'
gc_wig = 'summary_gc50.wig.gz'
vcf = 'test.vcf.gz'

ref_fasta = os.path.join(fasta_dir, ref_fasta)
bam_tumor = os.path.join(bam_dir, bam_tumor)
bam_normal = os.path.join(bam_dir, bam_normal)
pup_tumor = os.path.join(pup_dir, pup_tumor)
pup_normal = os.path.join(pup_dir, pup_normal)
gc_wig = os.path.join(gc_dir, gc_wig)
vcf = os.path.join(pup_dir, vcf)
# Some useful methods:


def fastq_bam(f1, f2, fasta, bam_out, head_id, bwa_bin='bwa',
              samtools_bin='samtools'):
    bwa_idx = "%s index %s" % (bwa_bin, fasta)
    bwa_mem = "%s mem -r '@RG\tID:%s\tSM:%s' %s %s %s" % (
        bwa_bin, head_id, head_id, fasta, f1, f2)
    bam_srt = "%s sort -l 9 -O bam --reference %s" % (
        samtools_bin, fasta)
    bam_idx = "%s index %s" % (
        samtools_bin, bam_out)
    bwa_idx = shlex.split(bwa_idx)
    bwa_mem = shlex.split(bwa_mem)
    bam_srt = shlex.split(bam_srt)
    bam_idx = shlex.split(bam_idx)
    proc1 = subprocess.Popen(bwa_idx)
    proc1.communicate()
    with open(bam_out, 'wb') as out:
        proc2 = subprocess.Popen(bwa_mem, stdout=subprocess.PIPE)
        proc3 = subprocess.Popen(bam_srt, stdin=proc2.stdout, stdout=out)
        proc3.communicate()
    proc4 = subprocess.Popen(bam_idx)
    proc4.communicate()


def fasta_idx(fasta, samtools_bin='samtools'):
    idx = "%s faidx %s" % (samtools_bin, fasta)
    idx = shlex.split(idx)
    proc1 = subprocess.Popen(idx)
    proc1.communicate()


def bgzip_file(file, bgzip_bin='bgzip', tabix_bin='tabix'):
    bgz = "%s %s" % (bgzip_bin, file)
    tbx = "%s -s 1 -b 2 -e 2 %s" % (tabix_bin, file + '.gz')
    bgz = shlex.split(bgz)
    tbx = shlex.split(tbx)
    proc1 = subprocess.Popen(bgz)
    proc1.communicate()
    proc2 = subprocess.Popen(tbx)
    proc2.communicate()


def vcf_mpileup(bam_normal, bam_tumor, output, fasta, samtools_bin='samtools'):
    samtools_mpileup = ('%s mpileup --reference %s '
                        '-v -t DP,AD -I -o %s %s %s') % (
        samtools_bin, fasta, output, bam_normal, bam_tumor)
    samtools_mpileup = shlex.split(samtools_mpileup)
    proc1 = subprocess.Popen(samtools_mpileup)
    proc1.communicate()

# Prepare the files


try:
    rmtree(tmp_dir)
    os.mkdir(tmp_dir)
except OSError:
    os.mkdir(tmp_dir)

# fasta index


fasta_idx(ref_fasta)

# BAM files
if not os.path.isdir(bam_dir):
    os.makedirs(bam_dir)
    fastq_bam(os.path.join(fastq_dir, 'testn_1.fq.gz'),
              os.path.join(fastq_dir, 'testn_2.fq.gz'),
              ref_fasta, bam_normal, 'sequenza_test_normal')
    fastq_bam(os.path.join(fastq_dir, 'testt_1.fq.gz'),
              os.path.join(fastq_dir, 'testt_2.fq.gz'),
              ref_fasta, bam_tumor, 'sequenza_test_tumor')


# mpileup files
if not os.path.isdir(pup_dir):
    os.makedirs(pup_dir)
    with xopen(pup_normal, 'wt') as normal_out:
        pileup_normal = bam_mpileup(bam_normal, ref_fasta)
        for line in pileup_normal:
            normal_out.write(line)
    bgzip_file(pup_normal)
    with xopen(pup_tumor, 'wt') as tumor_out:
        pileup_tumor = bam_mpileup(bam_tumor, ref_fasta)
        for line in pileup_tumor:
            tumor_out.write(line)
    bgzip_file(pup_tumor)

# VCF files
if not os.path.isfile(vcf):
    vcf_mpileup(bam_normal, bam_tumor, vcf, ref_fasta)

# gc_files

if not os.path.isdir(gc_dir):
    os.makedirs(gc_dir)
    with xopen(ref_fasta, 'rt') as fa_file, xopen(gc_wig, 'wt') as wig_out:
        stream_fasta = Fasta(fa_file, 50)
        wiggle = Wiggle(wig_out)
        for seq in stream_fasta:
            if seq:
                nucleotides = seq[3]
                Ns = countN(nucleotides, 'N')
                if Ns < 50:
                    if (seq[2] - seq[1] + 1) == 50:
                        gc = int(round(countN(nucleotides, 'G'), 0) +
                                 round(countN(nucleotides, 'C'), 0))
                        wiggle.write((seq[0], seq[1], (seq[2] + 1), gc))
