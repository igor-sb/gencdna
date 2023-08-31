"""CLI tool for SAM to alignment table."""

import fire
import pandas as pd

from gencdna.wgs.sam import (
    convert_genome_to_read_alignments,
    filter_gapped_exons,
    read_sam_to_dataframe,
)


def find_read_alignments(
    input_alignment_sam: str,
    input_read_positions_csv: str,
    output_alignment_table_csv: str,
    max_gap_len: int = 20_000,
) -> None:
    alignments_df = read_sam_to_dataframe(input_alignment_sam)
    reads_coords_df = pd.read_csv(input_read_positions_csv)
    alignments_reads_df = (
        alignments_df
        .pipe(convert_genome_to_read_alignments, reads_coords_df)
        .pipe(filter_gapped_exons, max_gap_len)
    )
    alignments_reads_df.to_csv(output_alignment_table_csv, index=False)


if __name__ == '__main__':
    fire.Fire(find_read_alignments)
