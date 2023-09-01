"""Test fixtures."""
import os
import pytest

fixture_path = '{base_path}/tests/exons_from_gff3/fixtures'.format(
    base_path=os.path.abspath('.'),
)

snapshot_path = '{base_path}/tests/exons_from_gff3/snapshots'.format(
    base_path=os.path.abspath('.'),
)


@pytest.fixture()
def sample_gff3():
    return '{fixture_path}/sample.gff3'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def ref_sample_exons_bed():
    return '{snapshot_path}/sample_exons.bed'.format(
        snapshot_path=snapshot_path,
    )


@pytest.fixture()
def ref_sample_annotations_csv():
    return '{snapshot_path}/sample_annotations.csv'.format(
        snapshot_path=snapshot_path,
    )
