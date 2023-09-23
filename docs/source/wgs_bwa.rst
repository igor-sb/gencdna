WGS BWA pipeline
================

The goal of this pipeline is to find gencDNAs among WGS reads. We look for
gencDNAs by aligning exons to reads then inspecting if the gap between two
subsequent exons is zero.

This pipeline proceeds by constructing a "genome" out of all reads in a sample,
then indexing that sample genome. With ``bwa index`` running single-threaded
(no idea how to multithread this?), this takes about 10 hours per sample.
However, multiple samples can be run in parallel.

I switched from Bowtie2 to BWA since the Bowtie2 pipeline did not work for all
genes. Some exons caused bowtie to search for the alignment for very long time
(I do not know why) and it never seemed to finish. BWA MEM does not have this
issue.

Circular consensus
------------------

Run circular consensus tool ``ccs``, where the default parameters are:

- ``--min-passes``: 3
- ``--min-snr``: 2.5
- ``--top-passes``: 60
- ``--min-rq``: 0.99 (min. predicted accuracy)

.. code-block:: bash

    ccs filename.bam filename.subreads.fastq

This is the only QC - because we want to try not to miss any rare reads that
may contain gencDNA.


Sample FASTQ to FASTA
---------------------

To create an index, we convert the FASTQ to a FASTA file:

.. code-block:: bash

    python gencdna/file_api/fastq_to_fasta.py \
        filename.subsread.fastq \
        reads.fasta


Sample genome index
-------------------

When supplied with a FASTA file containing multiple sequences, BWA does 
concatenation behind the scenes (and parses everything back), so we can build
an index directly from FASTA file form pre-processed reads:

.. code-block:: bash

    bwa index -p <bwa_index_prefix> reads.fasta


This is helpful, since I previously created code to do this manually with
bowtie2. 


Exon alignment
--------------


.. code-block:: bash

    bwa mem -t <num_cores> -a <bwa_index_prefix> <exons_fasta> > <exons_sam>


Filtering unmapped exons from SAM
---------------------------------

Unlike in bowtie2, there does not seem to be an option in ``bwa`` to omit
storing exons in the alignment SAM file that did not match any read. Therefore,
I added this step to filter those exons out:

.. code-block:: bash

    samtools view -F 4 -o <mapped_exons_sam> <exons_sam>


``-F`` flag filters out all reads that do not match 4 in the flag field. In bitwise
representation 4 represents unmapped reads. See page 7 at https://samtools.github.io/hts-specs/SAMv1.pdf
for details.


Alignment coordinates from SAM
------------------------------

Now that SAM file is free from exons that do not align to any read, we use
Python parse out SAM alignment and obtain alignment coordinates for each exon
on each read:

.. code-block:: bash

    python gencdna/file_api/alignment_coords.py \
        <mapped_exons_sam> \
        <mapped_exons_coords_csv>


Find exon joins
---------------

This script calculates gaps between adjacent exons then keeps the reads where
that gap is zero.

