"""Test fixtures."""

import os
import pytest

fixture_path = '{base_path}/tests/sam_parsing/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/sam_parsing/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def test_sam_file():
    return f'{fixture_path}/test.sam'


@pytest.fixture()
def ref_sam_parsed_to_csv_file():
    return f'{snapshot_path}/test.csv'
