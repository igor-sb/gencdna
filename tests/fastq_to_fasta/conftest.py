import os
import pytest


fixture_path = '{base_path}/tests/fastq_to_fasta/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/fastq_to_fasta/snapshots'.format(
    base_path=os.path.abspath('.'),
)



@pytest.fixture()
def reads_with_dups_fastq():
    return '{fixture_path}/reads_with_duplicates.fastq'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def expected_reads_with_dups_counts_csv():
    return '{snapshot_path}/reads_with_duplicates_counts.csv'.format(
        snapshot_path=snapshot_path,
    )


@pytest.fixture()
def expected_dedup_reads_fasta():
    return '{snapshot_path}/deduplicated_reads.fasta'.format(
        snapshot_path=snapshot_path,
    )


@pytest.fixture()
def reads_with_variable_lenghts_fasta():
    return '{fixture_path}/reads_with_variable_lengths.fasta'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def expected_reads_with_var_lenghts_counts_csv():
    return '{snapshot_path}/reads_with_var_lenghts_counts.csv'.format(
        snapshot_path=snapshot_path,
    )
