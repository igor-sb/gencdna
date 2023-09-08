"""Code for generating synthetic reads for alignment testing."""

from tests.setup.generate_sequences import generate_random_dna_sequence

LENGTH_FACTORS: dict[str, int] = {  # noqa: WPS407
    'SH': 10,
    'MD': 100,
    'LG': 1000,
    'VL': 10000,
}


def generate_random_exons(
    short_name: str,
    number_of_variants: int = 5,
) -> dict[str, str]:
    base_len = LENGTH_FACTORS[short_name]
    exons: dict[str, str] = {}
    for variant_id in range(1, number_of_variants + 1):
        for exon_id in range(1, 4):
            label = f'{short_name}{variant_id}E{exon_id}'
            exons[label] = generate_random_dna_sequence(
                base_len * exon_id,
            )
    return exons


def generate_random_genes(
    number_of_variants: int = 5,
) -> dict[str, str]:
    genes: dict[str, str] = {}
    for short_name in LENGTH_FACTORS.keys():
        genes.update(generate_random_exons(short_name, number_of_variants))
    return genes
