"""Test PCR primer filtering from a FASTA/FASTQ file."""

import os

from gencdna.file_api import filter_pcr_primers


def test_filter_reads_without_pcr_primers(
    mock_pcr_fasta_file,
    mock_pcr_filtered_fasta_file,
    snapshot,
):

    actual_fasta_output = filter_pcr_primers.filter_reads(
        input_file=mock_pcr_fasta_file,
        output_file='',
        forward_primer='ATGG',
        reverse_primer='GATT',
    )
    snapshot.snapshot_dir = os.path.dirname(
        mock_pcr_filtered_fasta_file,
    )
    snapshot.assert_match(
        actual_fasta_output,
        mock_pcr_filtered_fasta_file,
    )
