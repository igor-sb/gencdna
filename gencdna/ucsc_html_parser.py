"""Exon/intron parsing from URL."""

import io
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bs4 import BeautifulSoup


def parse_fasta_string(url_contents: str) -> str:
    soup = BeautifulSoup(url_contents, 'html.parser')
    pre_tag = soup.find('pre')
    if not pre_tag:
        raise ValueError('<pre> tag not found.')
    return pre_tag.text.strip()


def parse_fasta_records(
    url_contents: str,
) -> list[SeqRecord]:
    fasta_string = parse_fasta_string(url_contents)
    fasta_records: list[SeqRecord] = []
    with io.StringIO(fasta_string) as fasta:
        for fasta_record in SeqIO.parse(fasta, 'fasta'):
            fasta_records.append(fasta_record)
    return fasta_records


def parse_exons_introns(
    fasta_records: list[SeqRecord],
) -> tuple[list[SeqRecord], list[SeqRecord]]:
    exons, introns = [], []
    exon_index, intron_index = 0, 0
    for fasta_record in fasta_records:
        if re.match('^[ACGT]', str(fasta_record.seq)):
            exon_index += 1
            exons.append(SeqRecord(
                seq=fasta_record.seq,
                id='exon_{index:02}'.format(index=exon_index),
                name='',
                description='',
            ))
        elif re.match('^[acgt]', str(fasta_record.seq)):
            intron_index += 1
            introns.append(SeqRecord(
                seq=fasta_record.seq,
                id='intron_{index:02}'.format(index=intron_index),
                name='',
                description='',
            ))
    return (exons, introns)
