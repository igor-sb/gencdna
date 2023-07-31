"""Test blastn subprocess calls and output parsing."""
import os
import pytest

from pacbio_qc.blast import Blast, BlastOutputParser, is_executable_available


@pytest.mark.skipif(
    not is_executable_available('blastn'),
    reason='blastn not in path or not installed',
)
def test_blast_run(
    repeated_bases_flagged_fasta_file,
    target_with_repeated_bases_fasta_file,
    blast_output_from_repeated_bases_vs_target,
    snapshot,
):
    blast = Blast(
        query=repeated_bases_flagged_fasta_file,
        subject=target_with_repeated_bases_fasta_file,
    )
    blast.config['alignment']['word_size'] = 4
    blast_output = blast.run()
    blast_output_df = BlastOutputParser(blast_output).output_as_dataframe()
    snapshot.snapshot_dir = os.path.dirname(
        blast_output_from_repeated_bases_vs_target,
    )
    snapshot.assert_match(
        blast_output_df.to_csv(None, sep='\t', index=False),
        blast_output_from_repeated_bases_vs_target,
    )
