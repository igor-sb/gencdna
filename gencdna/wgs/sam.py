"""SAM alignment export to data frame."""

import numpy as np
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


def calculate_exon_gap(df):
    start_positions = df['read_start'].to_numpy()[1:]
    end_positions = df['read_end'].to_numpy()[:-1]
    start = np.concatenate(([np.nan], start_positions))
    end = np.concatenate(([np.nan], end_positions))
    df.loc[:, 'exon_gap'] = start - end
    return df


def filter_gapped_exons(
    alignments_reads_df: pd.DataFrame,
    max_gap_len: int,
) -> pd.DataFrame:
    alignments_with_gaps_df = (
        alignments_reads_df
        .groupby('read_id', group_keys=False)
        .apply(calculate_exon_gap)
        .sort_values(['read_id', 'read_start'])
    )
    exon_gaps = alignments_with_gaps_df['exon_gap']
    filter_mask = exon_gaps.isna() | (exon_gaps <= max_gap_len)  # noqa: WPS465
    return alignments_with_gaps_df[filter_mask]
