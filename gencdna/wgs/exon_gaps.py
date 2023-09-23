"""Find gaps between exons from the alignment coordinates table."""

import numpy as np
import pandas as pd


def parse_sequence_id_and_strand(exons_df: pd.DataFrame) -> pd.DataFrame:
    exons_df[['sequence_id', 'strand']] = (  # noqa: WPS359
        exons_df['query_id']
        .str.replace(')', '')
        .str.split('(', n=1, expand=True)
    )
    exons_df['sequence_id'] = exons_df['sequence_id'].astype(int)
    return exons_df


def rename_columns(exons_df: pd.DataFrame) -> pd.DataFrame:
    return (
        exons_df
        .rename(columns={
            'query_start': 'exon_start',
            'query_end': 'exon_end',
            'subject_id': 'read_id',
            'subject_start': 'read_start',
            'subject_end': 'read_end',
        })
        .drop(columns=['query_id'])
    )


def has_multiple_unique_values(df, column) -> bool:
    return df[column].nunique() > 1 and len(df) > 1


def find_reads_with_multiple_distinct_exons(
    exons_df: pd.DataFrame,
) -> pd.DataFrame:
    return (
        exons_df
        .groupby('read_id')
        .filter(lambda df: has_multiple_unique_values(df, 'read_start'))
    )


def calculate_exon_gaps(exons_df: pd.DataFrame) -> pd.DataFrame:
    # Calculate all gaps including between different read_ids for first exon
    exons_df['exon_gap'] = (
        exons_df['read_start'] - exons_df['read_end'].shift(1)
    )
    exons_df['index'] = exons_df.index
    first_indexes = exons_df.groupby('read_id').first()['index']

    # Correct the first gap to NaN
    exons_df['exon_gap'] = (
        exons_df['exon_gap']
        .where(~exons_df['index'].isin(first_indexes), np.NaN)
    )
    return exons_df.drop(columns=['index'])


def find_zero_gap_joins(exons_df: pd.DataFrame) -> pd.DataFrame:
    exons_pre_join = exons_df['exon_gap'].shift(-1) == 0
    exons_post_join = exons_df['exon_gap'] == 0
    return exons_df[exons_pre_join | exons_post_join]
