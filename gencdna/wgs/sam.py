"""SAM alignment export to data frame."""

import fire
import pandas as pd
from Bio.Align import sam

from gencdna.wgs.sam_records import record_to_dataframe


def read_sam_to_dataframe(sam_file: str) -> pd.DataFrame:
    aln_dfs = [
        record_to_dataframe(record)
        for record in sam.AlignmentIterator(sam_file)
    ]
    return pd.concat(aln_dfs)


def _genome_coord_to_reads_coord(alns_df, reads_df, which):
    return pd.merge(
        alns_df,
        reads_df[['read_id', which, 'read_interval']].rename(
            columns={
                'read_id': 'genome_{which}_read_id'.format(which=which),
                'start': 'read_start_on_genome',
                'read_interval': 'genome_{which}_interval'.format(which=which),
            },
        ),
        on='genome_{which}_interval'.format(which=which),
    )


def _remove_unused_columns(alns_reads_df: pd.DataFrame) -> pd.DataFrame:
    return alns_reads_df.drop(columns=[
        'genome_start_interval',
        'genome_end_interval',
        'genome_end_read_id',
        'end',
        'genome_id',
        'genome_start',
        'genome_end',
        'read_start_on_genome',
    ])


def _merge_genome_alignments_with_read_coords(
    alignments_df: pd.DataFrame,
    reads_coords_df: pd.DataFrame,
) -> pd.DataFrame:
    read_intervals_on_genome = pd.IntervalIndex.from_arrays(
        reads_coords_df['start'],
        reads_coords_df['end'],
        closed='both',
    )
    reads_coords_df['read_interval'] = read_intervals_on_genome
    alignments_df['genome_start_interval'] = pd.cut(
        alignments_df['genome_start'],
        read_intervals_on_genome,
    )
    alignments_df['genome_end_interval'] = pd.cut(
        alignments_df['genome_end'],
        read_intervals_on_genome,
    )
    return (
        alignments_df
        .pipe(_genome_coord_to_reads_coord, reads_coords_df, 'start')
        .pipe(_genome_coord_to_reads_coord, reads_coords_df, 'end')
        .query('genome_start_read_id == genome_end_read_id')
    )


def convert_genome_to_read_alignments(
    alignments_df: pd.DataFrame,
    reads_coords_df: pd.DataFrame,
) -> pd.DataFrame:

    alignments_with_reads_df = _merge_genome_alignments_with_read_coords(
        alignments_df,
        reads_coords_df,
    )

    return (
        alignments_with_reads_df
        .rename(columns={'genome_start_read_id': 'read_id'})
        .assign(
            read_start=lambda df: df.genome_start - df.read_start_on_genome,
            read_end=lambda df: df.genome_end - df.read_start_on_genome,
        )
        .sort_values(['read_id', 'read_start'])
        .pipe(_remove_unused_columns)
        .reset_index(drop=True)
    )


def find_read_alignments(
    input_alignment_sam: str,
    input_read_positions_csv: str,
    output_alignment_table_csv: str,
) -> None:
    alignments_df = read_sam_to_dataframe(input_alignment_sam)
    reads_coords_df = pd.read_csv(input_read_positions_csv)
    convert_genome_to_read_alignments(
        alignments_df,
        reads_coords_df,
    ).to_csv(output_alignment_table_csv, index=False)


if __name__ == '__main__':
    fire.Fire(find_read_alignments)
