Testing gencDNA finding using synthetic reads
=============================================

Implementation
--------------

In order to check the limitations of :ref:`wgs-bwa-pipeline`, and other tools
for detection insertion sites in the host genome, I created a code that
generates synthetic exon-exon joins of various lengths. 

The code is contained in ``gencdna/synthetic_reads/`` and is centered around the
class called ``ReadElement``. Each synthetic read is a obtained by concatenating
various read elements - which are either exons or introns.

You can use read elements like this in Python, where the first two arguments
are ``label`` and ``sequence``:

.. code-block:: python

        a = ReadElement('seq1', 'ATG')
        b = ReadElement('seq2', 'CCTA')
        print(a + b)
        #> seq1seq2 (7): ATGCCTA

so basically ``a + b`` is a new ``ReadElement`` object with label ``seq1seq2``
and sequence ``ATGCCTA``.


Usage
-----

To create a file for testing you can use:

.. code-block:: bash

    python gencdna/file_api/create_synthetic_reads.py \
        test_synthetic_reads.fasta \
        test_synthetic_exons.fasta

For reproducibility, the default seed is 1, but can be changed using
``--seed`` argument. This creates two files.

1. ``test_synthetic_exons.fasta`` file shows synthetic exon sequences.

::

    >[SH0E1:10]
    ATTCCCGTAA
    >[MD0E1:100]
    TCTACGATTAAGTCACAAC...
    >[LG0E1:1000]
    GCAAAGTTTTCTGAAATAA
    ...

2. ``test_synthetic_reads.fasta`` file shows synthetic read sequences.

::

    >-10000-[SH0E1:10][SH0E1:10]-10000-
    GTTTGGCATAATATACCGGAGGCTAGGCGTTATGTG...
    >-10000-[SH0E1:10][MD0E1:100]-100-
    CAGTTCAGAGTAAACGCCTGCAACCAATGCTAAGGA...
    >-10-[SH0E1:10][LG0E1:1000]-10000-
    CTACCGTAGGATTCCCGTAAGCAAAGTTTTCTGAA


Explanation of the output
-------------------------

The meta-data lines describe each read, where ``-length-`` represent introns of
length ``length``, and ``[name:length]`` represent exons named ``name`` and
having length ``length``. For example, ``-10000-[SH0E1:10][SH0E1:10]-1000-``
represents a read of total length 10000 + 10 + 10 + 1000 = 11020. Exons are
named as: ``{shortname}:{gene_variant}E{exon_number}``.

Genes and exons are generated from this hard-coded table of exon names and
lengths:

+------------+-----------+--------+--------+--------+--------+
| name       | shortname | exon_1 | exon_2 | exon_3 | intron |
+============+===========+========+========+========+========+
| short      | SH        | 10	  | 20     | 30     | 10     |
+------------+-----------+--------+--------+--------+--------+
| medium     | MD        | 100    | 200    | 300    | 100    |
+------------+-----------+--------+--------+--------+--------+
| long       | LG        | 1000   | 2000   | 3000   | 1000   |
+------------+-----------+--------+--------+--------+--------+
| verylong   | VL        | 10000  | 20000  | 30000  | 10000  |
+------------+-----------+--------+--------+--------+--------+

Basically there are 4 synthetic genes named SH, MD, LG and VL. Each has 3
exons of lengths that are multiple of 10 (short), 100 (medium), 1000 (long) or
10,000 (verylong). The sequence of each exon for each gene are generated
randomly, but reproducibly when given the same ``seed``. Introns are just
random sequences that we do not care about for this analysis.

So, the read labelled as ``-10000-[SH0E1:10][MD0E1:100]-100-`` contains a
concatenated sequence of length 20030:

1. intron of length 10000
2. exon of length 10 of short gene variant 0
3. exon of length 100 of medium gene variant 0
4. intron of length 100

And the sequences of all exons used in these reads are listed in the
``test_synthetic_reads.fasta``. So if we want to know what those exon sequences
are we can check that file and see that that the short gene variant 0, exon 1
``SH0E1`` has a sequence ``ATTCCCGTAA``.
