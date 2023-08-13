"""Code for aligning all exons to each of the reads."""

import fire
import pandas as pd

from gencdna.external_calls import BinaryExecWithYamlArgs
from gencdna.blast_parsing import BlastOutputParser


def align_exons_vs_reads(
    exons_fasta_file: str,
    reads_fasta_file: str,
    blast_config: str = 'config/blast.yml',
) -> pd.DataFrame:
    blast = BinaryExecWithYamlArgs('blastn', blast_config)
    blast_output = blast.run(
        '-lcase_masking',
        '-query',
        exons_fasta_file,
        '-subject',
        reads_fasta_file,
    )
    return (
        BlastOutputParser(blast_output)
        .output_as_dataframe()
        .sort_values(by=['subject_id', 'query_id'])
    )


def main(
    input_exons_fasta: str,
    input_reads_fasta: str,
    output_blast: str,
    blast_config: str = 'config/blast.yml',
) -> None:
    blast_output_df = align_exons_vs_reads(
        exons_fasta_file=input_exons_fasta,
        reads_fasta_file=input_reads_fasta,
        blast_config=blast_config,
    )
    blast_output_df.to_csv(output_blast, sep='\t', index=False)


if __name__ == '__main__':
    fire.Fire(main)
