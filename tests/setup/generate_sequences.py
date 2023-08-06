"""Random DNA sequence and quality score generation."""

import random


def generate_random_dna_sequence(length):
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choices(bases, k=length))


def generate_random_quality_characters(length, min_score=0, max_score=93):
    quality_scores = random.choices(range(min_score, max_score + 1), k=length)
    return ''.join(
        quality_score_to_character(qscore) for qscore in quality_scores
    )


def quality_character_to_score(quality_character: str, offset=33) -> int:
    return ord(quality_character) - offset


def quality_score_to_character(quality_score: int, offset=33) -> str:
    return chr(quality_score + offset)
