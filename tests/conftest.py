"""Test fixtures."""

import os
import pytest

fixture_path = '{base_path}/tests/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def snapshot_dir():
    return '{snapshot_path}'.format(
        snapshot_path=snapshot_path,
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


@pytest.fixture()
def blast_output_from_repeated_bases_vs_target():
    return '{snapshot_path}/blast_out_repeatedbases_vs_target.tsv'.format(
        snapshot_path=snapshot_path,
    )
