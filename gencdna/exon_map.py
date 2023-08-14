"""Create exon maps, e.g. "1 2' 3" for read with exons 1, 2rc and 3."""

import re
import pandas as pd


def parse_number_from_exon_id(exon_id: str) -> int:
    return int(re.sub('^exon_([0-9]+)$', '\\1', exon_id))


def determine_orientation(row: pd.Series) -> str:
    if row['subject_start'] <= row['subject_end']:        
        return '+'
    return '-'


def create_exon_label(row: pd.Series) -> str:
    if row['exon_orientation'] == '+':
        return str(row['exon_num'])
    return "{exon_num}'".format(exon_num=row['exon_num'])


def swap_subject_start_with_stop(blast_out: pd.DataFrame) -> pd.DataFrame:
    rows_with_rc_exons: list[bool] = blast_out['exon_orientation'] == '-'
    blast_out.loc[
        rows_with_rc_exons,
        ['subject_start', 'subject_end'],
    ] = blast_out.loc[
        rows_with_rc_exons,
        ['subject_end', 'subject_start'],
    ].values
    return blast_out


def create_exon_map(blast_out: pd.DataFrame) -> pd.DataFrame:
    blast_out['exon_num'] = blast_out['query_id'].apply(
        parse_number_from_exon_id,
    )
    blast_out['exon_orientation'] = blast_out.apply(
        determine_orientation,
        axis=1,
    )
    blast_out['exon_label'] = blast_out.apply(create_exon_label, axis=1)
    blast_out = swap_subject_start_with_stop(blast_out)    
    return (
        blast_out
        .sort_values(by=['subject_start'])
        .groupby('subject_id')['exon_label']
        .apply(lambda x: ' '.join(x), )
        .reset_index()
        .rename(columns={'exon_label': 'exon_map'})
    )
