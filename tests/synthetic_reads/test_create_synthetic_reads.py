"""Test file API for creating synthetic reads."""

import os
from tempfile import NamedTemporaryFile

from gencdna.file_api.create_synthetic_reads import create_synthetic_reads


def test_create_synthetic_reads(
    ref_synthetic_reads_fasta,
    snapshot,
):
    with NamedTemporaryFile('wb') as actual_file:
        create_synthetic_reads(
            actual_file.name,
            seed=1,
            number_of_length_variants=1,
            number_of_blocks_per_read=1,
            number_of_exons_in_block=2,
            selected_exon='exon_1',
        )
        snapshot.snapshot_dir = os.path.dirname(ref_synthetic_reads_fasta)
        with open(actual_file.name, 'rt') as actual_reads:
            snapshot.assert_match(
                actual_reads.read(),
                ref_synthetic_reads_fasta,
            )
