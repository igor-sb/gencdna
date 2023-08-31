Download exons for all human genes
==================================

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

We want to organize it in a bed file:

