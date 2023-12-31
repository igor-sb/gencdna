"""Test exons vs reads alignments."""

import os

from gencdna.file_api.align_exons_to_reads import align_exons_vs_single_read


def test_align_exons_vs_reads(
    exons_fasta_file,
    single_read_fasta_file,
    blast_config,
    alignment_output,
    snapshot,
):
    actual_df = align_exons_vs_single_read(
        exons_fasta_file=exons_fasta_file,
        read_fasta_file=single_read_fasta_file,
        config=blast_config,
    )
    snapshot.snapshot_dir = os.path.dirname(alignment_output)
    snapshot.assert_match(
        actual_df.to_csv(None, index=False),
        alignment_output,
    )
