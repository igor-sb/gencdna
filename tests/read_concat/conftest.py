"""Test fixtures."""
import os
import pytest

fixture_path = '{base_path}/tests/read_concat/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/read_concat/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def mock_pcr_fasta_file():
    return '{fixture_path}/reads.fasta'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def unconcatenated_fastq():
    return '{fixture_path}/reads.fastq'.format(fixture_path=fixture_path)


@pytest.fixture()
def ref_concatenated_fasta():
    return '{snapshot_path}/concatenated.fasta'.format(
        snapshot_path=snapshot_path,
    )


@pytest.fixture()
def ref_concatenated_positions():
    return '{snapshot_path}/concatenated_positions.csv'.format(
        snapshot_path=snapshot_path,
    )
