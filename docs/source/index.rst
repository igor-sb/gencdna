gencDNA finder for PacBio CCS DNA sequencing
============================================

This Python library contains code for detecting genomic cDNA, aka gencDNA
from the publication `Lee et al. Nature 2018 
<https://www.nature.com/articles/s41586-018-0718-6>`_.

Installation
------------

You need to have Python with `Poetry <https://python-poetry.org/>`_ installed.
See `Poetry Installation <https://python-poetry.org/docs/#installation>`_ 
instructions. Then run this:

.. code-block:: bash

   git clone https://github.com/igor-sb/gencdna.git
   cd gencdna
   make install

This will create an indentical Python virtual environment as what was used to
create this code and greatly minimize possible errors.

Python library structure
------------------------

`gencdna` is a Python library with various tools used to find gencDNAs.


Running the code
----------------

There are couple of ways to use the code. From Python for example:

.. code-block:: python

   from gencdna.file_api.expected_error_filter import filter_reads_with_low_expected_errors

   filter_reads_with_low_expected_errors(
      input_fastq_file,
      output_fastq_file,
      maximum_expected_errors = 0.01,
   )

... or from the command line:


.. code-block:: bash

   python gencdna/file_api/expected_error_filter.py \
      input.fastq \
      output.fastq \
      0.01



.. toctree::
   :maxdepth: 1
   :caption: Contents:

   workflow
   wgs_bowtie
   wgs_bwa
   exons
   all_genes_exons



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
