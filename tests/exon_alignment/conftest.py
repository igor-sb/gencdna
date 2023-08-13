"""Test fixtures."""

import os

import pytest

fixture_path = '{base_path}/tests/exon_alignment/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/exon_alignment/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def blast_config():
    return '{fixture_path}/blast_config.yml'.format(fixture_path=fixture_path)


@pytest.fixture()
def exons_fasta_file():
    return '{fixture_path}/exons_14nt.fasta'.format(fixture_path=fixture_path)


@pytest.fixture()
def single_read_fasta_file():
    return '{fixture_path}/read1_with_repeat_base_flags.fasta'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def reads_fasta_file():
    return '{fixture_path}/reads_with_repeat_base_flags.fasta'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def blast_output():
    return '{snapshot_path}/read1_with_repeat_base_flags_blastout.csv'.format(
        snapshot_path=snapshot_path,
    )
