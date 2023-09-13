"""CLI for creating synthetic reads with exon-exon joins."""

import random

import fire

from gencdna.synthetic_reads import (
    create_intron_exonblock_intron_read,
    random_exons,
)


def main(
    output_fasta_filename: str,
    seed: int = 1,
    number_of_exons: int = 5,
    number_of_read_replicates: int = 3,
) -> None:
    rng = random.Random(seed)
    exons = random_exons(rng, number_of_exons)
    reads = []
    for _ in range(number_of_read_replicates):
        reads.append(
            create_intron_exonblock_intron_read(
                exons,
                rng,
                number_of_exon_blocks=1,
                number_of_exons_per_block=1,
            ),
        )
    for _ in range(number_of_read_replicates):
        reads.append(
            create_intron_exonblock_intron_read(
                exons,
                rng,
                number_of_exon_blocks=2,
                number_of_exons_per_block=3,
            ),
        )
    with open(output_fasta_filename, 'wt') as output_fasta:
        for read in reads:
            output_fasta.write('>{label}\n{sequence}\n'.format(
                label=read['label'],
                sequence=read['sequence'],
            ))


if __name__ == '__main__':
    fire.Fire(main)
