import pytest

from pacbio_qc.blast import is_executable_available, Blast

@pytest.mark.skipif(
    not is_executable_available('blastn'),
    reason='blastn not in path or not installed'
)
def test_blast_run():
    blast = Blast()
    # test repeated_bases_flagged.fasta (query) vs targets_repeated_bases.fasta (subject)