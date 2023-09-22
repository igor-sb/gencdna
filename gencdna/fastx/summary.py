"""Count unique reads in a FASTQ file and output a FASTA with counts."""

import pandas as pd
from Bio import SeqIO

from gencdna.fastx.io import open_fastx_or_fastxgz


def count_unique_sequences(fastx_file: str, file_type: str) -> pd.DataFrame:
    sequence_counts: dict[str, int] = {}
    with open_fastx_or_fastxgz(fastx_file) as fastx_handle:
        for record in SeqIO.parse(fastx_handle, file_type):
            sequence = str(record.seq)
            sequence_counts[sequence] = sequence_counts.get(sequence, 0) + 1
        sequence_counts_df = pd.DataFrame(
            list(sequence_counts.items()),
            columns=['sequence', 'count'],
        )
    return sequence_counts_df.sort_values(
        by=['count', 'sequence'],
        ascending=[False, True],
    )


def dump_sequence_counts_to_fasta(sequence_counts_df: pd.DataFrame) -> str:
    fasta_dump = []
    for row in sequence_counts_df.itertuples():
        fasta_dump.append(
            '>read_{sequence_id:05};count_{count}\n{sequence}\n'.format(
                sequence_id=row.Index,
                count=row.count,
                sequence=row.sequence,
            ),
        )
    return ''.join(fasta_dump)


def summarize_read_lengths(fastx_file: str, file_type: str) -> pd.DataFrame:
    read_length_counts: dict[int, int] = {}
    with open_fastx_or_fastxgz(fastx_file, 'rt') as fastx_handle:
        for fastx_record in SeqIO.parse(fastx_handle, file_type):
            read_length = len(str(fastx_record.seq))
            read_length_counts[read_length] = (
                read_length_counts.get(read_length, 0) + 1
            )
    return pd.DataFrame({
        'read_length': read_length_counts.keys(),
        'count': read_length_counts.values(),
    }).sort_values(by=['count'], ascending=False)
