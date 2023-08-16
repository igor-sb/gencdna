"""Test fixtures."""

import os

import pytest

fixture_path = '{base_path}/tests/genome_browser/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/genome_browser/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def app_gene_html():
    return '{fixture_path}/app.html'.format(fixture_path=fixture_path)


@pytest.fixture()
def ref_app_gene_exons_and_introns_fasta():
    return '{snapshot_path}/app_exons_and_introns.fasta'.format(
        snapshot_path=snapshot_path,
    )
