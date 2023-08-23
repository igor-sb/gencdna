"""Generate perfect-quality FASTQ from FASTA (for bowtie2)."""

import fire
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def convert(
    input_fasta_file: str,
    output_fastq_file: str,
    max_qual: int = 93,
) -> None:

    reads: list[SeqRecord] = []
    for read in SeqIO.parse(input_fasta_file, 'fasta'):
        read.letter_annotations['phred_quality'] = [
            max_qual for _ in range(len(read.seq))
        ]
        reads.append(read)

    SeqIO.write(reads, output_fastq_file, 'fastq')


if __name__ == '__main__':
    fire.Fire(convert)
