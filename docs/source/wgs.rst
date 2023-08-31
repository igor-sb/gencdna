WGS pipeline
============

Preprocess
----------

Run ``ccs``:

.. code-block:: bash

    ccs filename.bam filename_ccs.fastq


Bowtie2 index
-------------

From each read create a bowtie2 index. We first concatenate all the reads into
a single FASTA file:

.. code-block:: bash

    python gencdna/file_api/wgs_reads.py \
        reads.fastq.gz \
        concat_reads.fasta \
        concat_positions.csv


Then build the index using the ``-a`` flag which is used to keep all the alignments:

.. code-block:: bash

    bowtie2-build \
        gencdna/concatenated.fasta \
        folder/concat_index


Alignments
----------

.. code-block:: bash

    bowtie2 -a --no-unal -L 10 \
        -x bowtie_index \
        -f exons.fasta \
        -S alignments.sam 


Tabulate results
----------------

.. code-block:: bash

    python gencdna/wgs/sam.py \
        alignments.sam \
