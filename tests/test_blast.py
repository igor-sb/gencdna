"""Test blastn subprocess calls and output parsing."""
import os

from gencdna.bash import BinaryExecWithYamlArgs
from gencdna.blast import BlastOutputParser


def test_blast_run(
    repeated_bases_flagged_fasta_file,
    target_with_repeated_bases_fasta_file,
    blast_output_from_repeated_bases_vs_target,
    snapshot,
):
    blast = BinaryExecWithYamlArgs('blastn', 'config/blast.yml')
    blast.config['arguments']['-word_size'] = 4
    blast_output = blast.run(
        '-query',
        repeated_bases_flagged_fasta_file,
        '-subject',
        target_with_repeated_bases_fasta_file,
    )
    blast_output_df = BlastOutputParser(blast_output).output_as_dataframe()
    snapshot.snapshot_dir = os.path.dirname(
        blast_output_from_repeated_bases_vs_target,
    )
    snapshot.assert_match(
        blast_output_df.to_csv(None, sep='\t', index=False),
        blast_output_from_repeated_bases_vs_target,
    )
