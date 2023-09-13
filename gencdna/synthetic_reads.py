"""Generate synthetic reads with exon-exon joins."""

import random
from typing import Optional

import pandas as pd

from tests.setup.random_seqs import random_sequence

# Used for generating mock exons. We have 4 genes, short, medium, long and
# verylong. Each will have exons of lengths:
#    base_length * 1, base_length * 2, base_length * 3, ...
EXON_LENGTHS: pd.DataFrame = pd.DataFrame({
    'name': ['short', 'medium', 'long', 'verylong'],
    'shortname': ['SH', 'MD', 'LG', 'VL'],
    'exon_1': [10, 100, 1000, 10000],
    'exon_2': [20, 200, 2000, 20000],
    'exon_3': [30, 300, 3000, 30000],
})

INTRON_LENGTHS: tuple[int, int, int, int] = (10, 100, 1000, 10000)


def create_intron_exonblock_intron_read(
    exons: list[tuple[str, str]],
    rng: random.Random,
    intron_length: Optional[int] = None,
    number_of_exon_blocks: int = 1,
    number_of_exons_per_block: int = 1,
) -> dict[str, str]:
    introns_and_exons = [random_intron(rng, intron_length)]

    for _ in range(number_of_exon_blocks):
        introns_and_exons.extend([
            *rng.choices(exons, k=number_of_exons_per_block),
            random_intron(rng, intron_length),
        ])
    read = {'label': '', 'sequence': ''}
    for intron_or_exon in introns_and_exons:
        read['label'] += intron_or_exon[0]
        read['sequence'] += intron_or_exon[1]
    return read


def random_intron(
    rng: random.Random,
    length: Optional[int] = None,
) -> tuple[str, str]:
    if not length:
        length = rng.sample(INTRON_LENGTHS, 1)[0]
    return (
        f'-{length}-',
        random_sequence(rng, length),
    )


def all_random_exons(rng: random.Random) -> list[tuple[str, str]]:
    genes_long = pd.melt(
        EXON_LENGTHS,
        id_vars=['name', 'shortname'],
        value_name='length',
        var_name='exon',
    )
    genes_long['exon_id'] = genes_long['exon'].str.replace('exon_', '')


def random_exons(
    rng: random.Random,
    number_of_exons: int = 1,
    shortname: Optional[str] = None,
) -> list[tuple[str, str]]:
    return [
        random_exon(rng, exon_number, shortname)
        for exon_number in range(number_of_exons)
    ]


def random_exon(
    rng: random.Random,
    gene_suffix: int,
    shortname: Optional[str] = None,
) -> tuple[str, str]:
    exon = select_random_exon(rng, shortname)
    return (
        '[{shortname}{number}E{exon_id}]'.format(
            shortname=exon['shortname'],
            number=gene_suffix,
            exon_id=exon['exon_id'],
        ),
        random_sequence(rng, length=int(exon['length'])),
    )


def select_random_exon(
    rng: random.Random,
    shortname: Optional[str] = None,
) -> dict[str, str]:
    genes_long = pd.melt(
        EXON_LENGTHS,
        id_vars=['name', 'shortname'],
        value_name='length',
        var_name='exon',
    )
    genes_long['exon_id'] = genes_long['exon'].str.replace('exon_', '')
    if shortname:
        genes_long = genes_long.query(f'shortname == "{shortname}"')
    random_row = rng.randrange(len(genes_long))
    return genes_long.iloc[random_row, :].to_dict()
