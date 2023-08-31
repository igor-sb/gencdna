"""Test phred score manipulation functions."""

import os
import tempfile

from gencdna.human_assembly import gff3_to_bed


def test_gff3_to_bed(
    annotation_gff3,
    ref_annotation_exons_bed,
    snapshot,
    bed_score=1000,
):
    with tempfile.NamedTemporaryFile() as actual_annotation_exons_bed_file:
        gff3_to_bed(
            input_gff3_file=annotation_gff3,
            output_bed_file=actual_annotation_exons_bed_file.name,
            bed_score=bed_score,
        )
        with open(actual_annotation_exons_bed_file.name, 'rt') as actual_bed:
            snapshot.snapshot_dir = os.path.dirname(ref_annotation_exons_bed)
            snapshot.assert_match(
                actual_bed.read(),
                ref_annotation_exons_bed,
            )
