"""CLI for FASTQ to FASTA aggregation with unique counts."""

import logging

import fire
import pandas as pd

from gencdna.count_unique_reads import (
    count_unique_sequences_in_fastq,
    dump_sequence_counts_to_fasta,
)

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


def write_unique_reads_from_fastq_to_fasta(
    input_fastq_file: str,
    output_fasta_file: str,
) -> None:
    sequence_counts: pd.DataFrame = count_unique_sequences_in_fastq(
        fastq_file=input_fastq_file,
    )
    with open(output_fasta_file, 'w') as output_fasta:
        output_fasta.write(dump_sequence_counts_to_fasta(sequence_counts))


if __name__ == '__main__':
    fire.Fire(write_unique_reads_from_fastq_to_fasta)
