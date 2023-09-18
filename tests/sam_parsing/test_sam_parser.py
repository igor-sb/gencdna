"""Test SAM coordinate parsing."""

import pandas as pd
from pandas.testing import assert_frame_equal

from gencdna.sam import alignment_coordinates, read_sam


def test_alignment_coordinates(test_sam_file, ref_sam_parsed_to_csv_file):
    sam_df = read_sam(test_sam_file)
    actual_df = alignment_coordinates(sam_df)
    ref_df = pd.read_csv(ref_sam_parsed_to_csv_file)
    assert_frame_equal(actual_df, ref_df)
