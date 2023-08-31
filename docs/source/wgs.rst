WGS pipeline
============

Commands and file names shown here are examples. Until I have put together
a full workflow use sample IDs are filenames or prefix to filenames.

Preprocess
----------

Run circular consensus tool ``ccs``, where the default parameters are:

- ``--min-passes``: 3
- ``--min-snr``: 2.5
- ``--top-passes``: 60
- ``--min-rq``: 0.99 (min. predicted accuracy)

.. code-block:: bash

    ccs filename.bam filename_ccs.fastq


Bowtie2 index
-------------

From each read create a bowtie2 index. We first concatenate all the reads into
a single FASTA file, ignoring the quality scores:

.. code-block:: bash

    python gencdna/wgs/concatenate.py \
        reads.fastq.gz \
        concatenated_reads.fasta \
        concatenated_positions.csv


Then build the index using the ``-a`` flag which is used to keep all the alignments:

.. code-block:: bash

    bowtie2-build \
        gencdna/concatenated_reads.fasta \
        folder/concatenated_reads


Align exon sequences to the bowtie2 index
-----------------------------------------

.. code-block:: bash

    bowtie2 -a --no-unal -L 10 \
        -x bowtie_index \
        -f exons.fasta \
        -S alignments.sam 

Alignment notes:

- ``-a`` keeps all alignments (by default on the best one is kept)
- ``--no-unal`` does not keep track of unaligned reads
- ``-f`` is the filename of the input fasta file with exons
- ``-S`` is the filename of the output SAM file with alignments


Tabulate results
----------------

This takes the SAM alignment file and ``concatenated_positions.csv`` as inputs
and creates an output ``alignment_table.csv`` that shows which exon was 
successfully aligned to which read.

.. code-block:: bash

    python gencdna/wgs/alignment_table.py \
        alignments.sam \
        concatenated_positions.csv \
        alignment_table.csv

The table will look like this:

.. csv-table:: 
    :file: alignment_table.csv
    :header-rows: 1


Filter results
--------------

At this point, to detect exon-exon joins, we look through the ``alignment_table.csv``
entries for each ``read_id``. If there are multiple exons in each ``read_id``,
and their ``exon_gap`` is low (ideally zero), then this would be an exon-exon join.

Keep in mind that the ``exon_gap`` is ``exon_start`` on that line minus
``exon_end`` on the previous, line. So in this case: 6745 - 309 = 6436.