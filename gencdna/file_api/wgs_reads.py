"""API for handling WGS reads."""

import re

import pandas as pd
from Bio import SeqIO

from gencdna.fastx_io import open_fastx_or_fastxgz


def concatenate(
    input_fastx: str,
    output_fasta: str,
    output_positions: str,
    output_annotation: str = 'concatenated_reads',
) -> None:
    positions: list[pd.DataFrame] = []
    start_position: int = 0
    with open(output_fasta, 'wt') as out_fa:
        out_fa.write('>{id}\n'.format(id=output_annotation))
        with open_fastx_or_fastxgz(input_fastx) as input_fx:
            for read in SeqIO.parse(input_fx, fastx_extension(input_fastx)):
                out_fa.write(str(read.seq))
                positions.append(pd.DataFrame({
                    'read_id': [read.id],
                    'start': [start_position],
                    'end': [start_position + len(read.seq) - 1],
                }))
                start_position += len(read.seq)
            out_fa.write('\n')
    pd.concat(positions).to_csv(output_positions, index=False)


def fastx_extension(filename: str) -> str:
    regex = re.search(r'.*\.(fast[a|q])(?:\.gz)?$', filename)
    if not regex:
        raise ValueError(
            'Invalid extension: {filename}'.format(filename=filename),
        )
    return regex.group(1)
