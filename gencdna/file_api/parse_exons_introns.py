"""CLI for exons/introns fasta from url."""

import fire
from Bio import SeqIO

from gencdna.ucsc.ucsc_html_parser import (
    parse_exons_introns,
    parse_fasta_records,
)


def write_exons_introns(
    ucsc_html: str,
    exons_fasta: str,
    introns_fasta: str = '',
    prefix: str = '',
) -> None:
    with open(ucsc_html) as url_contents:
        fasta_records = parse_fasta_records(url_contents.read())
        exons, introns = parse_exons_introns(fasta_records, prefix)
        SeqIO.write(exons, exons_fasta, 'fasta')
        if introns_fasta:
            SeqIO.write(introns, introns_fasta, 'fasta')


if __name__ == '__main__':
    fire.Fire(write_exons_introns)
