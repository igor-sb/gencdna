"""CLI for filtering exon alignments that are gapped / not exon-exon joins."""

import logging

import fire
import pandas as pd

from gencdna.wgs.exon_gaps import (
    calculate_exon_gaps,
    find_reads_with_multiple_distinct_exons,
    parse_sequence_id_and_strand,
    rename_columns,
)

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


def filter_gapped_alignments(
    mapped_exons_csv: str,
    exon_gaps_csv: str,
) -> None:
    LOG.info(f'Reading input: {mapped_exons_csv}')
    exons_df = pd.read_csv(mapped_exons_csv)
    LOG.info('Parsing sequence ID and strand info')
    exons_df = (
        exons_df
        .pipe(parse_sequence_id_and_strand)
        .pipe(rename_columns)
        .sort_values(['read_id', 'read_start'])
        .reset_index(drop=True)
    )
    LOG.info('Finding reads with multiple distinct exons')
    exons_df = (
        exons_df
        .pipe(find_reads_with_multiple_distinct_exons)
        .pipe(calculate_exon_gaps)
    )
    LOG.info('Saving output')
    exons_df.to_csv(exon_gaps_csv, index=False)


if __name__ == '__main__':
    fire.Fire(filter_gapped_alignments)
