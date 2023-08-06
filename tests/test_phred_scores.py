"""Test phred score manipulation functions."""

import gzip
import io
import json
import os

from Bio import SeqIO

from pacbio_qc.phred_scores import (
    expected_number_of_errors,
    filter_fastq_reads_by_expected_errors,
)


def test_expected_number_of_errors(example_fastq_file, snapshot):
    actual_expected_errors = []
    with gzip.open(example_fastq_file, 'rt') as fq:
        for fq_record in SeqIO.parse(fq, 'fastq'):
            actual_expected_errors.append(
                '{expected_error:.2E}'.format(
                    expected_error=expected_number_of_errors(fq_record),
                ),
            )
    snapshot.snapshot_dir = 'tests/snapshots'
    snapshot.assert_match(
        json.dumps(actual_expected_errors),
        'example_ccs_expected_errors.json',
    )


def test_filter_fastq_reads_by_expected_errors(
    example_reads_with_exons_fastq_file,
    example_reads_with_exons_filtered_fastq_file,
    snapshot,
):
    filtered_records = filter_fastq_reads_by_expected_errors(
        example_reads_with_exons_fastq_file,
        maximum_expected_errors=0.01,
    )
    with io.StringIO() as fastq_output:
        snapshot.snapshot_dir = os.path.dirname(
            example_reads_with_exons_filtered_fastq_file,
        )
        SeqIO.write(filtered_records, fastq_output, 'fastq')
        snapshot.assert_match(
            fastq_output.getvalue(),
            example_reads_with_exons_filtered_fastq_file,
        )
