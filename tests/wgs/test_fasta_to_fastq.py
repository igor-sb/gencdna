"""Test FASTA to FASTQ conversion."""

import os
import tempfile

from gencdna.wgs.fasta_to_fastq import convert


def test_fasta_to_fastq_convert(
    unconcatenated_fasta,
    ref_fastq_reads,
    snapshot,
    max_qual=93,
):
    with tempfile.NamedTemporaryFile('wt') as fastq_reads:
        convert(
            input_fasta_file=unconcatenated_fasta,
            output_fastq_file=fastq_reads.name,
            max_qual=max_qual,
        )
        snapshot.snapshot_dir = os.path.dirname(ref_fastq_reads)
        with open(fastq_reads.name, 'rt') as actual_fastq_reads:
            snapshot.assert_match(
                actual_fastq_reads.read(),
                ref_fastq_reads,
            )
