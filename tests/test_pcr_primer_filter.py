"""Test PCR primer filtering from a FASTA/FASTQ file."""

import os

from pacbio_qc.pcr_primer_filter import find_reads_with_pcr_primers


def test_filter_reads_without_pcr_primers(
    mock_pcr_experiment_fasta_file,
    mock_pcr_experiment_pcr_filtered_fasta_file,
    snapshot,
):

    actual_fasta_output = find_reads_with_pcr_primers(
        mock_pcr_experiment_fasta_file,
        cutadapt_args={
            '-g': 'ATGG...GATT',
            '--trimmed-only': '',
            '--minimum-length': '1',
            '--revcomp': '',
            '--quiet': '',
        },
    )
    snapshot.snapshot_dir = os.path.dirname(
        mock_pcr_experiment_pcr_filtered_fasta_file,
    )
    snapshot.assert_match(
        actual_fasta_output,
        mock_pcr_experiment_pcr_filtered_fasta_file,
    )
