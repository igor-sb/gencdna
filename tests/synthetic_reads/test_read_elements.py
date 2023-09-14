"""Tests for testing synthetic exons."""

import os
import pickle
import random
import re
from tempfile import NamedTemporaryFile

from gencdna.synthetic_reads.introns_and_exons import (
    create_synthetic_exons,
    random_synthetic_intron,
)


def test_create_synthetic_intron(intron_length=136):
    intron = random_synthetic_intron(random.Random(1), intron_length)
    sequence_dna_match = re.match('[ATGC]', intron.sequence)
    sequence_length_match = len(intron.sequence) == intron_length
    assert sequence_length_match and sequence_dna_match


def test_create_synthetic_exons(
    ref_synthetic_exons,
    snapshot,
):
    rng = random.Random(1)
    actual_exons = create_synthetic_exons(rng, 2)
    with NamedTemporaryFile('wb') as actual_exons_pkl:
        pickle.dump(actual_exons, actual_exons_pkl)
        snapshot.snapshot_dir = os.path.dirname(ref_synthetic_exons)
        with open(actual_exons_pkl.name, 'rb') as actual_pkl:
            snapshot.assert_match(
                actual_pkl.read(),
                ref_synthetic_exons,
            )
