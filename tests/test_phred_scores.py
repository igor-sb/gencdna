"""Test phred score manipulation functions."""

import gzip
import pickle

import pytest
from Bio import SeqIO

from pacbio_qc.phred_scores import expected_number_of_errors


@pytest.fixture(name='example_fastq_file')
def fixture_example_fastq_file():
    return 'tests/fixtures/example_ccs.fastq.gz'


def test_expected_number_of_errors(example_fastq_file, snapshot):
    actual_expected_errors = []
    with gzip.open(example_fastq_file, 'rt') as fq:
        for fq_record in SeqIO.parse(fq, 'fastq'):
            actual_expected_errors.append(
                expected_number_of_errors(fq_record),
            )
    snapshot.snapshot_dir = 'tests/snapshots'
    snapshot.assert_match(
        pickle.dumps(actual_expected_errors),
        'example_ccs_expected_errors.pkl',
    )
