"""Mock test URL download."""

from unittest.mock import Mock, patch
import pytest
from gencdna.ucsc.requests import exons_and_introns_fasta_from_url


@pytest.fixture(name='mock_requests_get')
def fixture_mock_requests_get():
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.content = b'<html><pre>FASTA content</pre></html>'
    return Mock(return_value=mock_response)


def test_exons_and_introns_fasta_from_url(mock_requests_get):
    with patch('requests.get', mock_requests_get):
        fasta_content = exons_and_introns_fasta_from_url('https://genome.ucsc.edu/cgi-bin/hgc')
        assert fasta_content == 'FASTA content'
