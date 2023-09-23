"""CLI for filtering FASTQ files by expected errors."""

import logging

import fire
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gencdna.fastx.expected_error_filter import (
    filter_fastq_reads_by_expected_errors,
)
from gencdna.fastx.io import open_fastx_or_fastxgz

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


def filter_reads_with_low_expected_errors(
    input_fastq_file: str,
    output_fastq_file: str,
    maximum_expected_errors: float = 1.0,
) -> None:
    filtered_records: list[SeqRecord] = filter_fastq_reads_by_expected_errors(
        input_fastq_file,
        maximum_expected_errors=maximum_expected_errors,
    )
    with open_fastx_or_fastxgz(output_fastq_file, 'wt') as output_fastq:
        SeqIO.write(filtered_records, output_fastq, 'fastq')


if __name__ == '__main__':
    fire.Fire(filter_reads_with_low_expected_errors)
