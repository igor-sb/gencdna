"""Test fixtures."""

import os
import pytest

fixture_path = '{base_path}/tests/repeat_bases/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/repeat_bases/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def repeated_bases_fastq_file():
    return '{fixture_path}/reads.fastq'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def repeated_bases_flagged_fastq_file():
    return '{snapshot_path}/flagged_reads.fastq'.format(
        snapshot_path=snapshot_path,
    )


@pytest.fixture()
def flagged_alignment_table_file():
    return '{fixture_path}/alignment_table.csv'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def ref_filled_alignment_table_file():
    return '{snapshot_path}/alignment_table_filled.csv'.format(
        snapshot_path=snapshot_path,
    )


@pytest.fixture()
def filled_alignment_table_file():
    return '{fixture_path}/alignment_table_filled.csv'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def ref_filled_reads_fasta_file():
    return '{snapshot_path}/flagged_reads_filled.fasta'.format(
        snapshot_path=snapshot_path,
    )
