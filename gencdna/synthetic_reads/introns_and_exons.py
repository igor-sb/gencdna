"""Code for creating synthetic exons and introns."""

import random
from typing import Optional

import pandas as pd

from gencdna.synthetic_reads.read_element import ReadElement
from tests.setup.random_seqs import random_sequence

# This table contains lengths of synthetic elements. We use these to
# programmatically generate random sequences of these lengths.
SYNTHETIC_ELEMENTS: pd.DataFrame = pd.DataFrame({
    'name': ['short', 'medium', 'long', 'verylong'],
    'shortname': ['SH', 'MD', 'LG', 'VL'],
    'exon_1': [10, 100, 1000, 10000],
    'exon_2': [20, 200, 2000, 20000],
    'exon_3': [30, 300, 3000, 30000],
    'intron': [10, 100, 1000, 10000],
})


def random_synthetic_intron(
    rng: random.Random,
    length: Optional[int] = None,
) -> ReadElement:
    length = length or random_synthetic_intron_length(rng)
    return ReadElement(
        f'-{length}-',
        random_sequence(rng, length),
    )


def random_synthetic_intron_length(
    rng: random.Random,
) -> str:
    return rng.sample(synthetic_intron_lengths(), 1)[0]


def create_synthetic_exons(
    rng: random.Random,
    number_of_length_variants: int = 1,
) -> list[ReadElement]:
    exons: list[ReadElement] = []
    for _, row in synthetic_exon_lengths().iterrows():
        for variant_number in range(number_of_length_variants):
            exons.append(
                ReadElement(
                    label=create_exon_label(row, variant_number),
                    sequence=random_sequence(rng, row.length),
                ),
            )
    return exons


def create_exon_label(
    row: pd.Series,
    variant_number: int,
) -> str:
    return '[{shortname}{variant}E{exon_id}:{length}]'.format(
        shortname=row['shortname'],
        variant=variant_number,
        exon_id=row['exon_id'],
        length=row['length'],
    )


def synthetic_exon_lengths() -> pd.DataFrame:
    exons_df = pd.melt(
        SYNTHETIC_ELEMENTS,
        id_vars=['name', 'shortname'],
        var_name='element',
        value_name='length',
    )
    return (
        exons_df
        .query('element != "intron"')
        .assign(
            exon_id=lambda df: df.element.str.replace('exon_', ''),
            length=lambda df: df.length.astype(int),
        )
        .sort_values('length')
    )


def synthetic_intron_lengths() -> pd.DataFrame:
    return SYNTHETIC_ELEMENTS['intron'].to_list()
