"""Test PCR primer filtering from a FASTA/FASTQ file."""

import os

from gencdna.file_api.pcr_primer_filter import (
    filter_reads_without_pcr_primers,
)


def test_filter_reads_without_pcr_primers(
    mock_pcr_experiment_fasta_file,
    mock_pcr_experiment_pcr_filtered_fasta_file,
    snapshot,
):

    actual_fasta_output = filter_reads_without_pcr_primers(
        input_file=mock_pcr_experiment_fasta_file,
        output_file='',
        forward_primer='ATGG',
        reverse_primer='GATT',
    )
    snapshot.snapshot_dir = os.path.dirname(
        mock_pcr_experiment_pcr_filtered_fasta_file,
    )
    snapshot.assert_match(
        actual_fasta_output,
        mock_pcr_experiment_pcr_filtered_fasta_file,
    )
