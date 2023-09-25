"""Parsers for human assembly FASTA/GFF3s."""

import pandas as pd
from gff3_parser import parse_gff3


def tabulate_exons_from_gff3(gff3_file: str) -> pd.DataFrame:
    return (
        parse_gff3(gff3_file, parse_attributes=True)
        .query('Type == "exon"')
        [[
            'Seqid',
            'Start',
            'End',
            'Strand',
            'gene_name',
            'transcript_name',
            'exon_number',
        ]]
        .rename(columns={'Seqid': 'chrom'})
    )


def add_annotations(df: pd.DataFrame) -> list[str]:
    df['Start'] = pd.to_numeric(df['Start'])
    df['End'] = pd.to_numeric(df['End'])
    df['exon_number'] = df['exon_number'].str.zfill(2)
    df['name'] = (
        df['transcript_name'] +  # noqa: W504
        ';' + df['gene_name'] +  # noqa: W504
        ';exon-' + df['exon_number']
    )
    df['sequence_id'] = (
        df
        .groupby(['chrom', 'Start', 'End'])
        .ngroup()
    )
    return df.sort_values(
        by=['chrom', 'Start', 'End', 'transcript_name'],
    )


def find_unique_exon_seqs_and_annotations(
    df: pd.DataFrame,
    bed_score: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = add_annotations(df)
    df = filter_short_and_long_exons(df)
    df['score'] = bed_score
    unique_seqs_df = (
        df[['chrom', 'Start', 'End', 'sequence_id', 'score', 'Strand']]
        .drop_duplicates()
    )
    annotation_df = (
        df[['sequence_id', 'gene_name', 'exon_number', 'transcript_name']]
        .drop_duplicates()
    )
    return (unique_seqs_df, annotation_df)


def filter_short_and_long_exons(
    df: pd.DataFrame,
    min_len: int = 30,
    max_len: int = 2000,
) -> pd.DataFrame:
    return (
        df
        .query('abs(End - Start) <= {max_len}'.format(max_len=max_len))
        .query('abs(End - Start) >= {min_len}'.format(min_len=min_len))
    )


def save_unique_exon_sequences_to_bed(
    unique_seqs_df: pd.DataFrame,
    output_bed_file: str,
) -> None:
    unique_seqs_df.to_csv(output_bed_file, sep='\t', header=None, index=False)


def save_unique_exon_annotations(
    annotation_df: pd.DataFrame,
    output_annotation_file: str,
) -> None:
    annotation_df.to_csv(output_annotation_file, sep='\t', index=False)
