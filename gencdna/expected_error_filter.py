"""Phred score manipulation functions."""

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from numpy.typing import NDArray

from gencdna.fastx_io import open_fastx_or_fastxgz


def convert_phred_scores_to_probability(phred_scores: NDArray) -> NDArray:
    return np.power(10, -np.array(phred_scores) / 10)


def expected_number_of_errors(fastq_record: SeqRecord) -> float:
    phred_scores = np.array(fastq_record.letter_annotations['phred_quality'])
    return np.sum(convert_phred_scores_to_probability(phred_scores))


def filter_fastq_reads_by_expected_errors(
    fastq_file: str,
    maximum_expected_errors: float = 1.0,
) -> list[SeqRecord]:
    filtered_fastq_records: list[SeqRecord] = []
    with open_fastx_or_fastxgz(fastq_file) as fastq_handle:
        for fastq_record in SeqIO.parse(fastq_handle, 'fastq'):
            number_of_errors = expected_number_of_errors(fastq_record)
            if number_of_errors <= maximum_expected_errors:
                filtered_fastq_records.append(fastq_record)
    return filtered_fastq_records
