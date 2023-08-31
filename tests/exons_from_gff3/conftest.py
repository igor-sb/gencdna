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
def annotation_gff3():
    return '{fixture_path}/annotation.gff3'.format(
        fixture_path=fixture_path,
    )


@pytest.fixture()
def ref_annotation_exons_bed():
    return '{snapshot_path}/annotation_exons.bed'.format(
        snapshot_path=snapshot_path,
    )
