"""HTTP retriever and parses for UCSC Genome Browser."""

import requests
from bs4 import BeautifulSoup


def exons_and_introns_fasta_from_url(
    url: str,
    timeout_seconds: int = 60,
) -> str:
    response = requests.get(url, timeout=timeout_seconds)
    if response.status_code == 200:  # noqa: WPS432
        soup = BeautifulSoup(response.content, 'html.parser')
        pre_tag = soup.find('pre')
        if pre_tag:
            fasta_content: str = pre_tag.text.strip()
        else:
            raise ValueError('<pre> tag not found.')
    else:
        raise IOError('Error: {code}'.format(code=response.status_code))
    return fasta_content
