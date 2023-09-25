"""CLI for flagging low quality repeated bases / homopolymer repeats."""

import logging

import fire
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gencdna.fastx.io import open_fastx_or_fastxgz
from gencdna.pcr.repeat_bases import (
    create_repeated_bases_regex,
    flag_lowqual_repeated_bases,
)

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


def flag_reads_with_low_quality_repeated_bases(
    input_fastq_file: str,
    output_fastq_file: str,
    quality_score_threshold: int = 30,
    flag_base: str = 'N',
    flag_quality_score: int = 93,
) -> None:
    output_records: list[SeqRecord] = []
    repeated_bases_regex = create_repeated_bases_regex()
    with open_fastx_or_fastxgz(input_fastq_file, 'rt') as input_fastq:
        for input_record in SeqIO.parse(input_fastq, 'fastq'):
            output_records.append(
                flag_lowqual_repeated_bases(
                    input_record,
                    repeated_bases_regex=repeated_bases_regex,
                    quality_score_threshold=quality_score_threshold,
                    flag_base=flag_base,
                    flag_quality_score=flag_quality_score,
                ),
            )
    with open_fastx_or_fastxgz(output_fastq_file, 'wt') as output_fastq:
        SeqIO.write(output_records, output_fastq, 'fastq')


if __name__ == '__main__':
    fire.Fire(flag_reads_with_low_quality_repeated_bases)
