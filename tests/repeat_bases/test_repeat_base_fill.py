"""Tests for filling in flagged repeat bases."""

import io
import os

from gencdna.file_api.fill_flagged_repeat_bases import (
    write_filled_alignment_table,
    write_filled_sequences_fasta,
)


def test_write_filled_alignment_table(
    flagged_alignment_table_file,
    ref_filled_alignment_table_file,
    snapshot,
):
    actual_alignment_table = io.StringIO()
    write_filled_alignment_table(
        flagged_alignment_table_file,
        actual_alignment_table,
        read_prefix='subject',
    )
    snapshot.snapshot_dir = os.path.dirname(
        ref_filled_alignment_table_file,
    )
    snapshot.assert_match(
        actual_alignment_table.getvalue(),
        ref_filled_alignment_table_file,
    )


def test_write_filled_sequences_fasta(
    filled_alignment_table_file,
    ref_filled_reads_fasta_file,
    snapshot,
):
    actual_fasta_file = io.StringIO()
    write_filled_sequences_fasta(
        filled_alignment_table_file,
        actual_fasta_file,
        read_prefix='subject',
    )
    snapshot.snapshot_dir = os.path.dirname(
        ref_filled_reads_fasta_file,
    )
    snapshot.assert_match(
        actual_fasta_file.getvalue(),
        ref_filled_reads_fasta_file,
    )
