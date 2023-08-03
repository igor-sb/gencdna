"""Quality control of homopolymer/repeated bases errors."""

import re
from typing import Any

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def flag_lowqual_repeated_bases(
    fastq_record: SeqRecord,
    repeated_bases_regex: re.Pattern,
    quality_score_threshold: int = 30,
    flag_base: str = 'N',
    flag_quality_score: int = 93,
) -> SeqRecord:
    """Flag low quality repeated bases (homopolymer repeats).

    Find repeated bases in a FASTQ record. If the first base has a quality
    score less than a threshold, flag the sequence by replacing the rest of the
    repeated bases with a flag_base character and assign it flag_quality_score.
    For example, replace ATTTTTG with ATNG, based on the quality score of the
    first T.
    """
    sequence = list(fastq_record.seq)
    quality_scores = fastq_record.letter_annotations['phred_quality']
    indexes_of_repeated_bases_after_second: list[int] = []

    for repeated_bases in repeated_bases_regex.finditer(str(fastq_record.seq)):
        first_quality_score_in_repeat = quality_scores[repeated_bases.start()]

        if first_quality_score_in_repeat < quality_score_threshold:
            sequence[repeated_bases.start() + 1] = flag_base
            quality_scores[repeated_bases.start() + 1] = flag_quality_score
            indexes_of_repeated_bases_after_second.extend(
                range(repeated_bases.start() + 2, repeated_bases.end()),
            )

    return construct_seq_record(
        metadata=fastq_record.id,
        sequence_list=remove_list_elements(
            sequence,
            indexes_of_repeated_bases_after_second,
        ),
        quality_scores_list=remove_list_elements(
            quality_scores,
            indexes_of_repeated_bases_after_second,
        ),
    )


def create_repeated_bases_regex(length: int = 2) -> re.Pattern:
    repeats = '{{{length},}}'.format(length=length)
    return re.compile('A{repeats}|C{repeats}|G{repeats}|T{repeats}'.format(
        repeats=repeats,
    ))


def remove_list_elements(
    input_list: list[Any],
    indices_to_remove: list[int],
) -> list[Any]:
    filtered_list = []
    for index, input_value in enumerate(input_list):
        if index not in indices_to_remove:
            filtered_list.append(input_value)
    return filtered_list


def construct_seq_record(metadata, sequence_list, quality_scores_list):
    return SeqRecord(
        seq=Seq(''.join(sequence_list)),
        id=metadata,
        description='',
        letter_annotations={'phred_quality': quality_scores_list},
    )
