"""Constructing synthetic reads for testing exon-exon join detection."""

import random

from tests.setup.generate_sequences import generate_random_dna_sequence


def generate_random_single_exon_reads(
    genes: dict[str, str],
    number_of_reads: int = 1,
    flank_lens: tuple[int, int, int, int] = (10, 100, 1000, 10000),
    seed: int = 42,
) -> dict[str, str]:
    random.seed(seed)
    reads: dict[str, str] = {}
    for read_id in range(number_of_reads):
        flanks = (
            generate_random_dna_sequence(random.sample(flank_lens, 1)[0]),
            generate_random_dna_sequence(random.sample(flank_lens, 1)[0]),
        )
        exon_label = random.sample(list(genes.keys()), 1)[0]
        read_label = make_read_label(
            read_id=read_id,
            left_flank=flanks[0],
            exon_label=exon_label,
            right_flank=flanks[1],
        )
        reads[read_label] = '{left_flank}{exon}{right_flank}'.format(
            left_flank=flanks[0],
            exon=genes[exon_label],
            right_flank=flanks[1],
        )
    return reads


def make_read_label(
    read_id: int,
    left_flank: str,
    exon_label: str,
    right_flank: str,
) -> str:
    left_flank_length = len(left_flank)
    right_flank_length = len(right_flank)
    read_structure = '-{left}-[{exon}]-{right}-'.format(
        left=left_flank_length,
        exon=exon_label,
        right=right_flank_length,
    )
    return f'read_{read_id}:{read_structure}'
