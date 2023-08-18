"""CLI for exons/introns fasta from url."""

import yaml
import fire
from Bio import SeqIO

from gencdna.genome_browser_api.parser import parse_exons_introns
from gencdna.genome_browser_api.url import (
    fasta_records_from_url,
    construct_url,
)


def write_exons_introns(
    gene_id: str,
    exons_fasta: str,
    introns_fasta: str,
    timeout_seconds: int = 60,
    geneid_file = 'config/genes.yml',
) -> None:
    url = construct_url(gene_id, config_file=geneid_file)
    fasta_records = fasta_records_from_url(url, timeout_seconds)
    exons, introns = parse_exons_introns(fasta_records)
    SeqIO.write(exons, exons_fasta, 'fasta')
    SeqIO.write(introns, introns_fasta, 'fasta')


if __name__ == '__main__':
    fire.Fire(write_exons_introns)
