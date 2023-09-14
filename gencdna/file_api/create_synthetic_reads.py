"""CLI for creating synthetic reads with exon-exon joins."""

import random

import fire

from gencdna.synthetic_reads.combinations import (
    create_intron_exonblocks_intron_elements,
)
from gencdna.synthetic_reads.introns_and_exons import create_synthetic_exons
from gencdna.synthetic_reads.read_element import ReadElement


def main(
    output_fasta_filename: str,
    seed: int = 1,
    number_of_exons: int = 5,
    number_of_read_replicates: int = 3,
) -> None:
    rng = random.Random(seed)
    exons = create_synthetic_exons(rng, number_of_exons)
    reads: list[ReadElement] = []
    for _ in range(number_of_read_replicates):
        reads.append(
            create_intron_exonblocks_intron_elements(
                exons,
                rng,
                number_of_blocks=1,
                number_of_exons_in_block=1,
            ),
        )
    for _ in range(number_of_read_replicates):
        reads.append(
            create_intron_exonblocks_intron_elements(
                exons,
                rng,
                number_of_blocks=2,
                number_of_exons_in_block=3,
            ),
        )
    with open(output_fasta_filename, 'wt') as output_fasta:
        for read in reads:
            output_fasta.write(read.to_fasta_str())


if __name__ == '__main__':
    fire.Fire(main)
