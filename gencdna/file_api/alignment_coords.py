"""CLI to extract alignment coordinates from SAM output."""

import fire

from gencdna.sam import alignment_coordinates, read_sam


def parse_alignment_coordinates(
    input_sam_file: str,
    output_csv_file: str,
) -> None:
    sam_df = read_sam(input_sam_file)
    align_coords_df = alignment_coordinates(sam_df)
    align_coords_df.to_csv(output_csv_file, index=False)


if __name__ == '__main__':
    fire.Fire(parse_alignment_coordinates)
