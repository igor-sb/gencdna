"""Generate synthetic reads with exon-exon joins."""

import random
from itertools import product
from typing import Any, Optional

from gencdna.synthetic_reads.introns_and_exons import random_synthetic_intron
from gencdna.synthetic_reads.read_element import (
    ReadElement,
    concat_read_elements,
)


def create_all_exonblocks_combinations(
    exons: list[ReadElement],
    number_of_exons_in_block: int = 1,
) -> list[ReadElement]:
    return [
        concat_read_elements(list(block_of_exons))
        for block_of_exons in product(exons, repeat=number_of_exons_in_block)
    ]


def partition_list(
    list_of_objects: list[Any],
    partition_size: int,
) -> list[list]:
    partitioned_list: list[list[Any]] = []
    for list_index in range(0, len(list_of_objects), partition_size):
        partitioned_list.append(
            list_of_objects[list_index:list_index + partition_size],
        )
    return partitioned_list


def concat_introns_to_ends(
    read_elements: list[ReadElement],
    rng: random.Random,
    intron_length: Optional[int] = None,
    ends: tuple[bool, bool] = (True, True),
) -> list[ReadElement]:
    elements: list[ReadElement] = []
    for element in read_elements:
        concat_element = ReadElement()
        if ends[0]:
            concat_element += random_synthetic_intron(rng, intron_length)
        concat_element += element
        if ends[1]:
            concat_element += random_synthetic_intron(rng, intron_length)
        elements.append(concat_element)
    return elements


def create_reads_with_exonblocks(
    exons: list[ReadElement],
    rng: random.Random,
    intron_length: Optional[int] = None,
    number_of_blocks_per_read: int = 1,
    number_of_exons_in_block: int = 1,
) -> list[ReadElement]:
    exon_blocks = create_all_exonblocks_combinations(
        exons,
        number_of_exons_in_block,
    )
    exon_intron_blocks = concat_introns_to_ends(
        exon_blocks,
        rng,
        intron_length,
        ends=(False, True),
    )
    exon_intron_blocks = partition_list(
        exon_intron_blocks,
        number_of_blocks_per_read,
    )
    reads = [concat_read_elements(blocks) for blocks in exon_intron_blocks]
    return concat_introns_to_ends(reads, rng, intron_length, (True, False))
