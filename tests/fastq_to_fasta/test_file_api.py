"""Test for counting unique reads in FASTQ or FASTA file."""

import os
import tempfile

from gencdna.file_api.fastq_unique_to_fasta import (
    write_unique_reads_from_fastq_to_fasta,
)


def test_write_unique_reads_from_fastq_to_fasta(
    reads_with_dups_fastq,
    expected_dedup_reads_fasta,
    snapshot,
):
    with tempfile.NamedTemporaryFile('wt') as actual_fasta_file:
        write_unique_reads_from_fastq_to_fasta(
            reads_with_dups_fastq,
            actual_fasta_file.name,
        )
        with open(actual_fasta_file.name, 'rt') as actual_fasta:
            snapshot.snapshot_dir = os.path.dirname(
                expected_dedup_reads_fasta,
            )
            snapshot.assert_match(
                actual_fasta.read(),
                expected_dedup_reads_fasta,
            )
