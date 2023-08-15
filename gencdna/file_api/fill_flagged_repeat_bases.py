"""CLI for filling in flagged repeat bases."""

import logging

import fire
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gencdna.repeat_bases import fill_alignment_output_sequences

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


def write_filled_alignment_output(
    flagged_alignment_table: str,
    filled_alignment_table: str,
    read_prefix: str = 'subject',
) -> None:
    align_out_df = pd.read_csv(flagged_alignment_table)
    filled_align_out_df = fill_alignment_output_sequences(
        align_out_df,
        read_prefix=read_prefix,
    )
    filled_align_out_df.to_csv(filled_alignment_table, index=False)


def write_filled_sequences_fasta(
    filled_alignment_table: str,
    fasta_file: str,
    read_prefix: str = 'subject',
) -> None:
    read_id = '{read_prefix}_id'.format(read_prefix=read_prefix)
    read_seq = '{read_prefix}_sequence'.format(read_prefix=read_prefix)
    filled_align_out_df = pd.read_csv(filled_alignment_table)
    fasta_records: list[SeqRecord] = filled_align_out_df.apply(
        lambda row: SeqRecord(Seq(row[read_seq]), row[read_id]),
    )
    SeqIO.write(fasta_records, fasta_file, 'fasta')


def write_outputs(
    input_alignment_table: str,
    output_alignment_table: str,
    output_fasta: str,
    read_prefix: str = 'subject',
) -> None:
    write_filled_alignment_output(
        input_alignment_table,
        output_alignment_table,
        read_prefix=read_prefix,
    )
    write_filled_sequences_fasta(
        output_alignment_table,
        output_fasta,
        read_prefix=read_prefix,
    )


if __name__ == '__main__':
    fire.Fire(write_outputs)
