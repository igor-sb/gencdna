"""Random DNA sequence and quality score generation."""

import random


def random_sequence(
    rng: random.Random,
    length: int,
):
    bases = ['A', 'C', 'G', 'T']
    return ''.join(rng.choices(bases, k=length))


def random_quality_characters(
    rng: random.Random,
    length: int,
    min_score: int = 0,
    max_score: int = 93,
):
    quality_scores = rng.choices(range(min_score, max_score + 1), k=length)
    return ''.join(
        quality_score_to_character(qscore) for qscore in quality_scores
    )


def quality_character_to_score(quality_character: str, offset=33) -> int:
    return ord(quality_character) - offset


def quality_score_to_character(quality_score: int, offset=33) -> str:
    return chr(quality_score + offset)
