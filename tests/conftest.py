"""Test fixtures."""

import os
import pytest

fixture_path = '{base_path}/tests/fixtures'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def example_fastq_file():
    return '{fixture_path}/example_ccs.fastq.gz'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def repeated_bases_fastq_file():
    return '{fixture_path}/repeated_bases.fastq'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def repeated_bases_flagged_fastq_file():
    return '{fixture_path}/repeated_bases_flagged.fastq'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def repeated_bases_flagged_fasta_file():
    return '{fixture_path}/repeated_bases_flagged.fasta'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def target_with_repeated_bases_fasta_file():
    return '{fixture_path}/targets_repeated_bases.fasta'.format(
        fixture_path=fixture_path,
    )
