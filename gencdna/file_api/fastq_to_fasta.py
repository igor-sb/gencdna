"""CLI for basic FASTQ to FASTA conversion."""

import fire
from Bio import SeqIO

from gencdna.fastx_io import open_fastx_or_fastxgz


def fastq_to_fasta(input_fastq_file: str, output_fasta_file: str) -> None:
    with open_fastx_or_fastxgz(input_fastq_file) as input_fastq:
        fastq_records = SeqIO.parse(input_fastq, 'fastq')
        SeqIO.write(fastq_records, output_fasta_file, 'fasta')


if __name__ == '__main__':
    fire.Fire(fastq_to_fasta)
