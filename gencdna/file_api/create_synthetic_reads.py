"""CLI for creating synthetic reads with exon-exon joins."""

import random
from typing import Optional

import fire

from gencdna.synthetic_reads.combinations import create_reads_with_exonblocks
from gencdna.synthetic_reads.introns_and_exons import create_synthetic_exons


def create_synthetic_reads(  # noqa: WPS211, WPS210
    output_reads_fasta: str,
    output_exons_fasta: str,
    seed: int = 1,
    number_of_length_variants: int = 1,
    number_of_blocks_per_read: int = 1,
    number_of_exons_in_block: int = 2,
    selected_exon: Optional[str] = None,
) -> None:
    rng = random.Random(seed)
    exons = create_synthetic_exons(
        rng,
        number_of_length_variants,
        selected_exon,
    )
    reads = create_reads_with_exonblocks(
        exons,
        rng,
        number_of_blocks_per_read=number_of_blocks_per_read,
        number_of_exons_in_block=number_of_exons_in_block,
    )
    with open(output_reads_fasta, 'wt') as reads_fasta:
        for read in reads:
            reads_fasta.write(read.to_fasta_str())

    with open(output_exons_fasta, 'wt') as exons_fasta:
        for exon in exons:
            exons_fasta.write(exon.to_fasta_str())


if __name__ == '__main__':
    fire.Fire(create_synthetic_reads)
