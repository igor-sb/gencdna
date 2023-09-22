"""CLI to extract alignment coordinates from SAM output."""

import logging

import fire

from gencdna.sam.sam import alignment_coordinates, read_sam

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


def parse_alignment_coordinates(
    input_sam_file: str,
    output_csv_file: str,
) -> None:
    LOG.info(f'Reading {input_sam_file}')
    sam_df = read_sam(input_sam_file)
    n_records = len(sam_df)
    LOG.info(f' {n_records} total records.')
    LOG.info('Calculating coordinates:')
    align_coords_df = alignment_coordinates(sam_df, LOG)
    LOG.info('Writing output')
    align_coords_df.to_csv(output_csv_file, index=False)


if __name__ == '__main__':
    fire.Fire(parse_alignment_coordinates)
