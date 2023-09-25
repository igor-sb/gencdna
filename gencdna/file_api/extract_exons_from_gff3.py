"""CLI for GFF3 to BED conversion and exon extraction."""

import fire

from gencdna.gff3.exon_parser import (
    find_unique_exon_seqs_and_annotations,
    tabulate_exons_from_gff3,
)


def gff3_to_bed(
    input_full_gff3: str,
    output_unique_seqs_bed: str,
    output_annotation_csv: str,
    bed_score: int = 1000,
) -> None:
    exons_df = tabulate_exons_from_gff3(input_full_gff3)
    unique_seqs_df, annotations_df = find_unique_exon_seqs_and_annotations(
        exons_df,
        bed_score,
    )
    unique_seqs_df.to_csv(
        output_unique_seqs_bed,
        sep='\t',
        header=None,
        index=False,
    )
    annotations_df.to_csv(output_annotation_csv, index=False)


if __name__ == '__main__':
    fire.Fire(gff3_to_bed)
