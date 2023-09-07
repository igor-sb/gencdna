Download and extract exons for all human genes
==============================================

First, download the Human Genome FASTA file and GFF3 file with annotations. We
will use Gencode release 44 from: https://www.gencodegenes.org/human/release_44.html
the two files are listed under:

- Basic gene annotation (chr): It contains the basic gene annotation on the reference chromosomes only (GFF3)
- Genome sequence, primary assembly (GRCh38) (FASTA)

.. code-block:: bash

    wget 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gff3.gz'
    wget 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz'
    gunzip gencode.v44.basic.annotation.gff3.gz
    gunzip GRCh38.primary_assembly.genome.fa.gz

The plan here is to convert the goofy GFF3 format, where gene names are buried in a string such as::

    ID=ENSG00000223972.6;gene_id=ENSG00000223972.6;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;level=2;hgnc_id=HGNC:37102;havana_gene=OTTHUMG00000000961.2

into a tab-delimited file or a csv where we can link the gene names, transcript
name, and other meta-data to the unique sequence on the chromosome. 

First, we will organize the chromosomal positions in a bed file and other meta-data
in a csv file:

.. code-block:: bash

    python gencdna/file_api/extract_exons_from_gff3.py \
        gencode.v44.basic.annotation.gff3 \
        all_exons.bed \
        all_exons_positions.csv

BED file is a simple tab-delimited file, without a header row, which contains
the columns:

- chromosome ID
- chromosome start position
- chromosome end position
- sequence ID (used to cross-reference FASTA and BED files)
- sequence score (used for UCSC Genome Browser plotting)
- strand (+ or -)

Then use `bedtools <https://bedtools.readthedocs.io/en/latest/>`_:

.. code-block:: bash

    bedtools getfasta -nameOnly -s \
        -fi GRCh38.p14.genome.fa \
        -bed all_exons.bed \
        -fo all_exons.fasta

This will use genomic locations from bed file and DNA sequences from a FASTA
file to pull out sequences in those genomic locations and store them in 
``all_exons.fasta``. This is the file we will use for bowtie2 alignment.

Resulting FASTA file shoud look like this::

    >0(+)
    TTAACTTGCCGT...
    >1(+)
    TGTCTGACTTCG...
    >2(-)
    TGGAGGAAAGAT...

Where (+) or (-) indicate the strand bedtools used and the number is the
sequence ID number in the BED file and in the ``all_exons_positions.csv`` file::

    sequence_idgene_name,exon_number,transcript_name
    0,DDX11L2,01,DDX11L2-202
    1,DDX11L1,01,DDX11L1-201
    2,DDX11L1,02,DDX11L1-201
    ...
