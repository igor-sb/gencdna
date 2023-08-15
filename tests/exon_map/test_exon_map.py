"""Test creating exon maps for reads."""

import pandas as pd
from pandas.testing import assert_frame_equal

from gencdna.exon_map import create_exon_map


def test_create_exon_map(alignment_output):
    alignment_output_df = pd.read_csv(alignment_output)
    actual_exon_map = (
        create_exon_map(alignment_output_df)
        .sort_values(by=['subject_id'])
    )
    reference_exon_map = pd.DataFrame({
        'subject_id': [
            'read_exons1,2,3_but_revcomp_exon2',
            'read_exons1,2,3_with_repeat_base_error',
            'read_exons1,2,3_with_repeat_base_flags_and_flanks',
        ],
        'exon_map': [
            "1 2' 3",
            '1 2 3',
            '1 2 3',
        ],
    })
    assert_frame_equal(actual_exon_map, reference_exon_map)
