"""Test fixtures."""

import os
import pytest

snapshot_path = '{base_path}/tests/synthetic_reads/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def ref_synthetic_exons():
    return f'{snapshot_path}/synthetic_exons.txt'


@pytest.fixture()
def ref_synthetic_reads_fasta():
    return f'{snapshot_path}/synthetic_reads.fasta'
