"""API for handling WGS reads."""

import pandas as pd
from Bio import SeqIO


def concatenate(
    input_fastq: str,
    output_fasta: str,
    output_positions: str,
    output_annotation: str = 'concatenated_reads',
) -> None:
    positions: list[pd.DataFrame] = []
    start_position: int = 0
    read_length: int = 0
    with open(output_fasta, 'wt') as out_fa:
        out_fa.write('>{id}\n'.format(id=output_annotation))
        for read in SeqIO.parse(input_fastq, 'fastq'):
            read_length = len(read.seq)
            out_fa.write(str(read.seq))
            positions.append(pd.DataFrame({
                'read_id': [read.id],
                'start': [start_position],
                'end': [start_position + read_length - 1],
            }))
            start_position += read_length
        out_fa.write('\n')
    pd.concat(positions).to_csv(output_positions, index=False)
