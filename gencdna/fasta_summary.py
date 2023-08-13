"""FASTA file statistics."""

import pandas as pd
from Bio import SeqIO

from gencdna.fastx_io import open_fastx_or_fastxgz


def summarize_read_lengths(fasta_file: str) -> pd.DataFrame:
    read_length_counts: dict[int, int] = {}
    with open_fastx_or_fastxgz(fasta_file, 'rt') as fasta_handle:
        for fasta_record in SeqIO.parse(fasta_handle, 'fasta'):
            read_length = len(str(fasta_record.seq))
            read_length_counts[read_length] = (
                read_length_counts.get(read_length, 0) + 1
            )
    return pd.DataFrame({
        'read_length': read_length_counts.keys(),
        'count': read_length_counts.values(),
    }).sort_values(by=['count'], ascending=False)
