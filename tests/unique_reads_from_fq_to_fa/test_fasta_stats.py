"""Tests for read length statistics in a FASTA file."""

import os

from gencdna.fasta_summary import summarize_fasta_read_lengths


def test_summarize_fasta_read_lengths(
    mock_pcr_experiment_fasta_file,
    mock_pcr_experiment_read_lengths_file,
    snapshot,
):
    actual_df = summarize_fasta_read_lengths(mock_pcr_experiment_fasta_file)
    snapshot.snapshot_dir = os.path.dirname(
        mock_pcr_experiment_read_lengths_file,
    )
    snapshot.assert_match(
        actual_df.to_csv(None, index=False),
        mock_pcr_experiment_read_lengths_file,
    )
