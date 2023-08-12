gencDNA finder workflow
=======================

1. Filter reads by the expected number of errors:

.. code-block:: bash

    poetry run python \
        gencdna/file_api/expected_error_filter.py \
        raw_reads.fastq.gz \
        filtered_reads.fastq.gz \
        0.01

2. Flag reads with low quality repeated bases (homopolymer errors):

.. code-block:: bash

    poetry run python

X. Search for each exon across reads then output a BLAST outfmt 6 output:

.. code-block:: bash

    poetry run python