"""Test phred score manipulation functions."""

import os
import re
import tempfile

import pytest
from Bio import SeqIO

from gencdna.expected_error_filter import expected_number_of_errors
from gencdna.file_api.filter_by_expected_error import (
    filter_reads_with_low_expected_errors,
)


def test_expected_number_of_errors(
    reads_with_known_errors,
    tolerance=0.00001,
):
    with open(reads_with_known_errors, 'rt') as fq:
        for fq_record in SeqIO.parse(fq, 'fastq'):
            actual_error = expected_number_of_errors(fq_record)
            record_id = str(fq_record.id)
            ref_error = float(
                re.sub(
                    '.*expected_error_(.*)$',
                    '\\1',  # noqa: WPS342
                    record_id,
                ),
            )
            assert ref_error == pytest.approx(actual_error, abs=tolerance)


def test_filter_fastq_reads_by_expected_errors(
    reads_with_known_errors,
    reads_with_known_errors_filtered,
    snapshot,
    max_errors=0.01,
):
    with tempfile.NamedTemporaryFile() as actual_filtered_reads_file:
        filter_reads_with_low_expected_errors(
            input_fastq_file=reads_with_known_errors,
            output_fastq_file=actual_filtered_reads_file.name,
            maximum_expected_errors=max_errors,
        )
        with open(actual_filtered_reads_file.name, 'rt') as actual_reads:
            snapshot.snapshot_dir = os.path.dirname(
                reads_with_known_errors_filtered,
            )
            snapshot.assert_match(
                actual_reads.read(),
                reads_with_known_errors_filtered,
            )
