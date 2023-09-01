"""Test phred score manipulation functions."""

import os
import tempfile

from gencdna.file_api.extract_exons_from_gff3 import gff3_to_bed


def test_gff3_to_bed(
    sample_gff3,
    ref_sample_exons_bed,
    ref_sample_annotations_csv,
    snapshot,
    bed_score=1000,
):
    with tempfile.NamedTemporaryFile() as actual_bed:
        with tempfile.NamedTemporaryFile() as actual_csv:
            gff3_to_bed(
                input_full_gff3=sample_gff3,
                output_unique_seqs_bed=actual_bed.name,
                output_annotation_csv=actual_csv.name,
                bed_score=bed_score,
            )
            with open(actual_bed.name, 'rt') as actual_bed_file:
                snapshot.snapshot_dir = os.path.dirname(ref_sample_exons_bed)
                snapshot.assert_match(
                    actual_bed_file.read(),
                    ref_sample_exons_bed,
                )
            with open(actual_csv.name, 'rt') as actual_csv_file:
                snapshot.snapshot_dir = os.path.dirname(
                    ref_sample_annotations_csv,
                )
                snapshot.assert_match(
                    actual_csv_file.read(),
                    ref_sample_annotations_csv,
                )
