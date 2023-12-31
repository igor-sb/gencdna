"""Test homopolymer/repeated bases quality control code."""

import os

from Bio import SeqIO

from gencdna.pcr.repeat_bases import (
    create_repeated_bases_regex,
    flag_lowqual_repeated_bases,
)


def test_flag_lowqual_repeated_bases(
    repeated_bases_fastq_file,
    repeated_bases_flagged_fastq_file,
    snapshot,
):
    output_fastq_strings = []
    regex = create_repeated_bases_regex()
    for input_fq_record in SeqIO.parse(repeated_bases_fastq_file, 'fastq'):
        output_fq_record = flag_lowqual_repeated_bases(input_fq_record, regex)
        output_fastq_strings.append(output_fq_record.format('fastq'))
    snapshot.snapshot_dir = os.path.dirname(
        repeated_bases_flagged_fastq_file,
    )
    snapshot.assert_match(
        ''.join(output_fastq_strings),
        repeated_bases_flagged_fastq_file,
    )
