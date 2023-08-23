"""Test read concatenation for WGS data."""

import os
import tempfile

from gencdna.wgs.concatenate import concatenate


def test_concate_fastq_reads(
    unconcatenated_fastq,
    ref_concatenated_fasta,
    ref_concatenated_positions,
    snapshot,
):
    with tempfile.NamedTemporaryFile('wt') as concatenated_fasta:
        with tempfile.NamedTemporaryFile('wt') as concatenated_positions:
            concatenate(
                input_fastx=unconcatenated_fastq,
                output_fasta=concatenated_fasta.name,
                output_positions=concatenated_positions.name,
                output_annotation='concatenated_read',
            )
            snapshot.snapshot_dir = os.path.dirname(
                ref_concatenated_fasta,
            )

            with open(concatenated_fasta.name, 'rt') as concatenated_fa:
                snapshot.assert_match(
                    concatenated_fa.read(),
                    ref_concatenated_fasta,
                )

            with open(concatenated_positions.name, 'rt') as concatenated_pos:
                snapshot.assert_match(
                    concatenated_pos.read(),
                    ref_concatenated_positions,
                )


def test_concate_fasta_reads(
    unconcatenated_fasta,
    ref_concatenated_fasta,
    ref_concatenated_positions,
    snapshot,
):
    with tempfile.NamedTemporaryFile('wt') as concatenated_fasta:
        with tempfile.NamedTemporaryFile('wt') as concatenated_positions:
            concatenate(
                input_fastx=unconcatenated_fasta,
                output_fasta=concatenated_fasta.name,
                output_positions=concatenated_positions.name,
                output_annotation='concatenated_read',
            )
            snapshot.snapshot_dir = os.path.dirname(
                ref_concatenated_fasta,
            )

            with open(concatenated_fasta.name, 'rt') as concatenated_fa:
                snapshot.assert_match(
                    concatenated_fa.read(),
                    ref_concatenated_fasta,
                )

            with open(concatenated_positions.name, 'rt') as concatenated_pos:
                snapshot.assert_match(
                    concatenated_pos.read(),
                    ref_concatenated_positions,
                )


def test_concate_fastq_gz_reads(
    unconcatenated_fastq_gz,
    ref_concatenated_fasta,
    ref_concatenated_positions,
    snapshot,
):
    with tempfile.NamedTemporaryFile('wt') as concatenated_fasta:
        with tempfile.NamedTemporaryFile('wt') as concatenated_positions:
            concatenate(
                input_fastx=unconcatenated_fastq_gz,
                output_fasta=concatenated_fasta.name,
                output_positions=concatenated_positions.name,
                output_annotation='concatenated_read',
            )
            snapshot.snapshot_dir = os.path.dirname(
                ref_concatenated_fasta,
            )

            with open(concatenated_fasta.name, 'rt') as concatenated_fa:
                snapshot.assert_match(
                    concatenated_fa.read(),
                    ref_concatenated_fasta,
                )

            with open(concatenated_positions.name, 'rt') as concatenated_pos:
                snapshot.assert_match(
                    concatenated_pos.read(),
                    ref_concatenated_positions,
                )
