"""CLI for GFF3 to BED conversion and exon extraction."""

import fire

from gencdna.exon_parser import (
    save_tabulated_exons_to_bed,
    tabulate_exons_from_gff3,
)


def gff3_to_bed(
    input_gff3_file: str,
    output_bed_file: str,
    bed_score: int = 1000,
) -> None:
    exons_df = tabulate_exons_from_gff3(input_gff3_file)
    save_tabulated_exons_to_bed(exons_df, output_bed_file, bed_score)


if __name__ == '__main__':
    fire.Fire(gff3_to_bed)
