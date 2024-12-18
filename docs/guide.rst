User cookbook
=============

In order to process BAM files and generate file index, the software `samtools`_ and `tabix`_
need to be installed in the system.

The package *sequenza-utils* includes several programs and it should support the generation
of *seqz* files using commonly available input files, such as fasta, BAM and vcf files.

To write your own program using the *sequenza-utils* library, lease refer to the
:doc:`API library interface <api>`


Generate GC reference file
--------------------------

The GC content source required to generate *seqz* files must to be in the  `wiggle track`_
format (WIG). In order to generate the wig file from any fasta file use the :code:`gc_wiggle`
program.

.. code:: bash

    sequenza-utils gc_wiggle --fasta genome.fa.gz -w 50 -o genome_gc50.wig.gz


From BAM files
--------------

Normal and tumor BAM files
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    sequenza-utils bam2seqz --normal normal_sample.bam --tumor tumor_sample.bam \
        --fasta genome.fa.gz -gc genome_gc50.wig.gz --output sample.seqz.gz


Normal and tumor pileup files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    sequenza-utils bam2seqz --normal normal_sample.pielup.gz \
        --tumor tumor_sample.pielup.gz --fasta genome.fa.gz \
        -gc genome_gc50.wig.gz --output sample.seqz.gz --pileup

Without normal, workaround
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    sequenza-utils bam2seqz --normal tumor_sample.bam --tumor tumor_sample.bam \
        --normal2 non_matching_normal_sample.bam --fasta genome.fa.gz \
        -gc genome_gc50.wig.gz --output sample.seqz.gz


Binning seqz, reduce memory
---------------------------


.. code:: bash

    sequenza-utils seqz_binning --seqz sample.seqz.gz --window 50 \
        -o sample_bin50.seqz.gz


Seqz from VCF files
-------------------

VCF files with **DP** and **AD** tags
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    sequenza-utils snp2seqz --vcf sample_calls.vcf.gz -gc genome_gc50.wig.gz \
        --output samples.seqz.gz


Mutect/Caveman/Strelka2 preset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    sequenza-utils snp2seqz --vcf sample_calls.vcf.gz -gc genome_gc50.wig.gz \
        --preset mutect --output samples.seqz.gz



.. code:: bash

    sequenza-utils snp2seqz --vcf sample_calls.vcf.gz -gc genome_gc50.wig.gz \
        --preset caveman --output samples.seqz.gz



.. code:: bash

    sequenza-utils snp2seqz --vcf sample_calls.vcf.gz -gc genome_gc50.wig.gz \
        --preset strelka2_som --output samples.seqz.gz




Merge seqz files
----------------


Non overlapping calls (eg different chromosomes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    gzcat sample_chr1.seqz.gz sample_chr1.seqz.gz | \
        gawk '{if (NR!=1 && $1 != "chromosome") {print $0}}' | bgzip > \
        sample.seqz.gz
    tabix -f -s 1 -b 2 -e 2 -S 1 sample.seqz.gz


Overlapping sample_calls
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    sequenza-utils seqz_merge --seqz1 sample_somatic.seqz.gz \
        --seqz2 sample_snps.seqz.gz --output samples.seqz.gz



.. _`samtools`: http://samtools.sourceforge.net
.. _`tabix`: http://www.htslib.org/doc/tabix.html
.. _`wiggle track`: https://genome.ucsc.edu/goldenpath/help/wiggle.html