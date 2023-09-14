"""Test functions for generating read element combinations."""

import random

from gencdna.synthetic_reads.combinations import (
    create_all_exonblocks_combinations,
    create_reads_with_exonblocks,
)
from gencdna.synthetic_reads.read_element import ReadElement


def test_create_all_exonblocks_combinations():
    exons = [
        ReadElement('exon1', 'ATG'),
        ReadElement('exon2', 'CCTA'),
    ]
    actual_combinations = create_all_exonblocks_combinations(exons, 2)
    ref_combinations = [
        ReadElement('exon1exon1', 'ATGATG'),
        ReadElement('exon1exon2', 'ATGCCTA'),
        ReadElement('exon2exon1', 'CCTAATG'),
        ReadElement('exon2exon2', 'CCTACCTA'),
    ]
    assert actual_combinations == ref_combinations


def test_create_intron_exonblocks_intron_elements():
    rng = random.Random(1)
    exons = [ReadElement('[exon1]', 'ATG'), ReadElement('[exon2]', 'CCTA')]
    actual_reads = create_reads_with_exonblocks(
        exons,
        rng,
        intron_length=5,
        number_of_blocks_per_read=2,
        number_of_exons_in_block=3,
    )
    ref_reads = [
        ReadElement(
            '-5-[exon1][exon1][exon1]-5-[exon1][exon1][exon2]-5-',
            'TTACGATGATGATGATTCCATGATGCCTACGTAA',
        ),
        ReadElement(
            '-5-[exon1][exon2][exon1]-5-[exon1][exon2][exon2]-5-',
            'GTCTGATGCCTAATGTCTACATGCCTACCTAGATTA',
        ),
        ReadElement(
            '-5-[exon2][exon1][exon1]-5-[exon2][exon1][exon2]-5-',
            'CGTTGCCTAATGATGAGTCACCTAATGCCTACAACC',
        ),
        ReadElement(
            '-5-[exon2][exon2][exon1]-5-[exon2][exon2][exon2]-5-',
            'GAATCCCTACCTAATGAAACCCCTACCTACCTAATGGA',
        ),
    ]
    assert actual_reads == ref_reads
