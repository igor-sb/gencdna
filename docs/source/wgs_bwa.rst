WGS BWA pipeline
================

The WGS Bowtie2 pipeline does not work for all genes. Some exons cause bowtie
to search for the alignment for very long time (I am not sure exactly why), but
BWA MEM does not have that issue and runs much quicker anyway.

Preprocess
----------

Run circular consensus tool ``ccs``, where the default parameters are:

- ``--min-passes``: 3
- ``--min-snr``: 2.5
- ``--top-passes``: 60
- ``--min-rq``: 0.99 (min. predicted accuracy)

.. code-block:: bash

    ccs filename.bam filename_ccs.fastq

BWA index
---------

When supplied with a FASTA file containing multiple sequences, BWA does 
concatenation behind the scenes (and parses everything back), so we can build
an index directly from FASTA file form pre-processed reads:

.. code-block:: bash

    bwa index -t <num_cores> -p <bwa_index_prefix> <reads_fasta>


BWA alignment
-------------


.. code-block:: bash

    bwa mem -t <num_cores> -a <bwa_index_prefix> <exons_fasta> > <output_sam>


Filter unmapped exons from SAM
------------------------------

.. code-block:: bash

    samtools view -F 4 -o <output_mapped_sam> <output_sam>


-F flag filters out all reads that do not match 4 in the flag field. In bitwise
representation 4 represents unmapped reads. See page 7 at https://samtools.github.io/hts-specs/SAMv1.pdf
for details.


Alignment coordinates from SAM
------------------------------

.. code-block:: bash

