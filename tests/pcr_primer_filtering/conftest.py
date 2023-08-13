import os
import pytest


fixture_path = '{base_path}/tests/pcr_primer_filtering/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/pcr_primer_filtering/snapshots'.format(
    base_path=os.path.abspath('.'),
)



@pytest.fixture()
def mock_pcr_fasta_file():
    return '{fixture_path}/reads.fasta'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def mock_pcr_filtered_fasta_file():
    return '{snapshot_path}/reads_filtered.fasta'.format(
        snapshot_path=snapshot_path,
    )