"""Mock test URL download."""

import io
import os
from unittest.mock import Mock, patch

import pytest
from Bio import SeqIO

from gencdna.genome_browser_api.url import fasta_records_from_url


@pytest.fixture(name='local_requests_get')
def fixture_local_requests_get(app_gene_html):
    mock_response = Mock()
    mock_response.status_code = 200
    with open(app_gene_html) as html_file:
        mock_response.content = html_file.read()
    return Mock(return_value=mock_response)


def test_exons_and_introns_fasta_from_url(
    local_requests_get,
    ref_app_gene_exons_and_introns_fasta,
    snapshot,
):
    with patch('requests.get', local_requests_get):
        fasta_records = fasta_records_from_url('https://genome.ucsc.edu/')
        with io.StringIO() as app_gene_exons_and_introns_fasta:
            snapshot.snapshot_dir = os.path.dirname(
                ref_app_gene_exons_and_introns_fasta,
            )
            SeqIO.write(
                fasta_records,
                app_gene_exons_and_introns_fasta,
                'fasta',
            )
            snapshot.assert_match(
                app_gene_exons_and_introns_fasta.getvalue(),
                ref_app_gene_exons_and_introns_fasta,
            )
