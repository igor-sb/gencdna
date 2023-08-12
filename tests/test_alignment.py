"""Test exons vs reads alignments."""

import os

from gencdna.file_api.exon_alignment import align_exons_vs_reads


def test_align_exons_vs_reads(
    example_exons_fasta_file,
    example_reads_with_exons_flagged_fasta_file,
    example_reads_with_exons_yml,
    example_reads_with_exons_blastout,
    snapshot,
):
    actual_df = align_exons_vs_reads(
        exons_fasta_file=example_exons_fasta_file,
        reads_fasta_file=example_reads_with_exons_flagged_fasta_file,
        blast_config=example_reads_with_exons_yml,
    )
    snapshot.snapshot_dir = os.path.dirname(
        example_reads_with_exons_blastout,
    )
    snapshot.assert_match(
        actual_df.to_csv(None, index=False),
        example_reads_with_exons_blastout,
    )
