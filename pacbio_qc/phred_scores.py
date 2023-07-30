"""Phred score manipulation functions."""

import Bio
import numpy as np
from numpy.typing import NDArray


def convert_phred_scores_to_probability(phred_scores: NDArray) -> NDArray:
    return np.power(10, -np.array(phred_scores) / 10)


def expected_number_of_errors(fastq_record) -> float:
    phred_scores = np.array(fastq_record.letter_annotations['phred_quality'])
    return np.sum(convert_phred_scores_to_probability(phred_scores))
