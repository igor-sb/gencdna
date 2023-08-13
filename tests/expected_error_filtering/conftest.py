"""Test fixtures."""
import os
import pytest

fixture_path = '{base_path}/tests/expected_error_filtering/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/expected_error_filtering/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def reads_with_known_errors():
    return '{fixture_path}/reads.fastq'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def reads_with_known_errors_filtered():
    return '{snapshot_path}/reads_filtered.fastq'.format(
        snapshot_path=snapshot_path,
    )
