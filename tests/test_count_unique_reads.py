"""Test for counting unique reads in FASTQ or FASTA file."""

import gzip
import os
import pandas as pd

from pacbio_qc.count_unique_reads import (
    count_unique_sequences_in_fastq,
    dump_sequence_counts_to_fasta,
)


def test_count_unique_sequences_in_fastq(
    example_fastq_file,
    example_ccs_sequence_counts_file,
    snapshot,
):
    with gzip.open(example_fastq_file, 'rt') as fq:
        actual_df = count_unique_sequences_in_fastq(fq)
    snapshot.snapshot_dir = os.path.dirname(
        example_ccs_sequence_counts_file,
    )
    snapshot.assert_match(
        actual_df.to_csv(None, index=False),
        example_ccs_sequence_counts_file,
    )


def test_dump_sequence_counts_to_fasta(
    example_ccs_sequence_counts_file,
    example_ccs_unique_reads_fasta_file,
    snapshot,
):
    sequence_counts_df = pd.read_csv(example_ccs_sequence_counts_file)
    actual_fasta = dump_sequence_counts_to_fasta(sequence_counts_df)
    snapshot.snapshot_dir = os.path.dirname(
        example_ccs_unique_reads_fasta_file,
    )
    snapshot.assert_match(
        actual_fasta,
        example_ccs_unique_reads_fasta_file,
    )
