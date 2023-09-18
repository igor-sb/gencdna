"""SAM file parsing into a table with 0-based coordinates."""

import re

import pandas as pd

OPERATIONS = {  # noqa: WPS407
    'M': 'matches_or_mismatches',
    'I': 'insertions',
    'D': 'deletions',
    'N': 'skips',
    'S': 'softclips',
    'H': 'hardclips',
    'P': 'pads',
    '=': 'seq_matches',
    'X': 'seq_mismatches',
}


class SamAlignment(object):

    def __init__(self, sam_record: pd.Series) -> None:
        self.cigars: list[str] = re.findall(
            r'\d+[MIDNSHP=X]',
            sam_record['cigar'],
        )
        self.reference_start_pos: int = int(sam_record['reference_start'])
        self.query_name: str = sam_record['query_name']
        self.reference_name: str = sam_record['reference_name']
        self.opers: dict[str, int] = {
            'matches_or_mismatches': 0,
            'insertions': 0,
            'deletions': 0,
            'skips': 0,
            'softclips': 0,
            'hardclips': 0,
        }

        for cigar in self.cigars:
            oper_type = OPERATIONS[cigar[-1]]
            self.opers[oper_type] = (
                self.opers.get(oper_type, 0) + int(cigar[:-1])
            )

    def read_length(self) -> int:
        return sum((
            self.opers['matches_or_mismatches'],
            self.opers['insertions'],
            self.opers['skips'],
            self.opers['softclips'],
            self.opers['hardclips'],
        ))

    def read_start(self) -> int:
        first_cigar = self.cigars[0]
        if first_cigar.endswith('S') or first_cigar.endswith('H'):
            return int(first_cigar[:-1])
        return 0

    def read_end(self) -> int:
        read_length: int = self.read_length()
        end_clip: int = 0
        last_cigar = self.cigars[-1]
        if last_cigar.endswith('S') or last_cigar.endswith('H'):
            end_clip = int(last_cigar[:-1])
        return read_length - end_clip - 1

    def reference_start(self) -> int:
        return self.reference_start_pos - 1

    def reference_end(self) -> int:
        return sum((
            self.reference_start(),
            self.opers['matches_or_mismatches'],
            self.opers['deletions'],
            self.opers['skips'],
        )) - 1

    def as_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame({
            'query_id': [self.query_name],
            'subject_id': [self.reference_name],
            'query_start': [self.read_start()],
            'query_end': [self.read_end()],
            'subject_start': [self.reference_start()],
            'subject_end': [self.reference_end()],
        })


def read_sam(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(
        filepath,
        comment='@',
        sep='\t',
        header=None,
        usecols=range(11),  # noqa: WPS432
        names=[
            'query_name',
            'bitwise_flags',
            'reference_name',
            'reference_start',
            'mapping_quality',
            'cigar',
            'reference_name_mate',
            'position_mate',
            'template_length',
            'query_sequence',
            'query_quality',
        ],
    )
    df['bitwise_flags'] = [f'{flag:012b}' for flag in df['bitwise_flags']]
    return df


def alignment_coordinates(sam_records: pd.DataFrame) -> pd.DataFrame:
    sam_records_coords: list[pd.DataFrame] = []
    for _, row in sam_records.iterrows():
        sam_records_coords.append(SamAlignment(row).as_dataframe())
    return pd.concat(sam_records_coords, ignore_index=True)
