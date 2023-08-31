"""Parsers for human assembly FASTA/GFF3s."""

import pandas as pd
from gff3_parser import parse_gff3


def tabulate_exons_from_gff3(gff3_file: str) -> pd.DataFrame:
    df = (
        parse_gff3(gff3_file, parse_attributes=True)
        .query('Type == "exon"')
    )
    return df[[
        'Seqid',
        'Start',
        'End',
        'Strand',
        'gene_name',
        'exon_number',
        'gene_type',
    ]]


def build_gene_exon_annotation(df: pd.DataFrame) -> list[str]:
    gene_exon_annotations: list[str] = []
    for gene, exon_num in zip(df['gene_name'], df['exon_number']):
        gene_exon_annotations.append(
            '{gene}_exon_{exon_num:02}'.format(
                gene=gene,
                exon_num=int(exon_num),
            ),
        )
    return gene_exon_annotations


def save_tabulated_exons_to_bed(
    df: pd.DataFrame,
    bed_file: str,
    bed_score: int = 1000,
) -> None:
    df['name'] = build_gene_exon_annotation(df)
    df['score'] = bed_score
    df = (
        df[['Seqid', 'Start', 'End', 'name', 'score', 'Strand']]
        .rename(
            columns={
                'Seqid': 'chrom',
                'Start': 'chromStart',
                'End': 'chromEnd',
                'Strand': 'strand',
            },
        )
    )
    df.to_csv(bed_file, sep='\t', header=False, index=False)
