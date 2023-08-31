"""Functions for parsing SAM records using Bio.Align."""

import pandas as pd
from Bio.Align import Alignment


def exon_coordinates(record: Alignment) -> tuple[int, int]:
    exon_positions = record.coordinates[1]
    return (exon_positions[0], exon_positions[-1])


def genome_coordinates(record: Alignment) -> tuple[int, int]:
    genome_positions = record.coordinates[0]
    return (genome_positions[0], genome_positions[-1])


def exon_start_end_revcomp(record: Alignment) -> tuple[int, int, bool]:
    coordinates = exon_coordinates(record)
    if coordinates[0] > coordinates[1]:
        return (coordinates[1], coordinates[0], True)
    return (coordinates[0], coordinates[1], False)


def record_to_dataframe(record: Alignment) -> pd.DataFrame:
    exon_start, exon_end, exon_rc = exon_start_end_revcomp(record)
    genome_start, genome_end = genome_coordinates(record)
    return pd.DataFrame({
        'exon_id': [record.query.id],
        'genome_id': [record.target.id],
        'exon_start': [exon_start],
        'exon_end': [exon_end],
        'exon_rc': [exon_rc],
        'genome_start': [genome_start],
        'genome_end': [genome_end],
        'exon_seq': [str(record.query.seq[exon_start:exon_end])],
        'genome_seq': [str(record.target.seq[genome_start:genome_end])],
    })
