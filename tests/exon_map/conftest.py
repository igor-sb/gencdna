"""Test fixtures."""

import os
import pytest

fixture_path = '{base_path}/tests/exon_map/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/exon_map/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def blast_output():
    return '{fixture_path}/blast_output.csv'.format(
        fixture_path=fixture_path,
    )
