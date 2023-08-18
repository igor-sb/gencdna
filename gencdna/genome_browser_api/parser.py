"""Exon/intron parsing from URL."""

import re

from Bio.SeqRecord import SeqRecord


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
