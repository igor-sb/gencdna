"""Test fixtures."""

import os
import pytest

fixture_path = '{base_path}/tests/repeated_base_flagging/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/repeated_base_flagging/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def repeated_bases_fastq_file():
    return '{fixture_path}/reads.fastq'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def repeated_bases_flagged_fastq_file():
    return '{snapshot_path}/reads_flagged.fastq'.format(
        snapshot_path=snapshot_path,
    )
