"""Tests for read length statistics in a FASTA file."""

import os

import pandas as pd

from gencdna.fastx.summary import (
    count_unique_sequences,
    dump_sequence_counts_to_fasta,
    summarize_read_lengths,
)


def test_count_unique_sequences_in_fastq(
    reads_with_dups_fastq,
    expected_reads_with_dups_counts_csv,
    snapshot,
):
    actual_df = count_unique_sequences(reads_with_dups_fastq, 'fastq')
    snapshot.snapshot_dir = os.path.dirname(
        expected_reads_with_dups_counts_csv,
    )
    snapshot.assert_match(
        actual_df.to_csv(None, index=False),
        expected_reads_with_dups_counts_csv,
    )


def test_dump_sequence_counts_to_fasta(
    expected_reads_with_dups_counts_csv,
    expected_dedup_reads_fasta,
    snapshot,
):
    sequence_counts_df = pd.read_csv(expected_reads_with_dups_counts_csv)
    actual_reads_fasta = dump_sequence_counts_to_fasta(sequence_counts_df)
    snapshot.snapshot_dir = os.path.dirname(
        expected_dedup_reads_fasta,
    )
    snapshot.assert_match(
        actual_reads_fasta,
        expected_dedup_reads_fasta,
    )


def test_summarize_read_lengths(
    reads_with_variable_lenghts_fasta,
    expected_reads_with_var_lenghts_counts_csv,
    snapshot,
):
    actual_df = summarize_read_lengths(
        reads_with_variable_lenghts_fasta,
        file_type='fasta',
    )
    snapshot.snapshot_dir = os.path.dirname(
        expected_reads_with_var_lenghts_counts_csv,
    )
    snapshot.assert_match(
        actual_df.to_csv(None, index=False),
        expected_reads_with_var_lenghts_counts_csv,
    )
