"""URL handling for Genome Browser HTTP GET API."""

import io
from typing import Any

import requests
import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bs4 import BeautifulSoup


def config_dict_to_get_string(config: dict[str, Any]) -> str:
    key_vals = [
        '{key}={val}'.format(key=get_key, val=get_val)
        for get_key, get_val in zip(config.keys(), config.values())
    ]
    return '&'.join(key_vals)


def construct_url(
    gene_id: str,
    chromosome: str,
    left_position: int,
    right_position: int,
    config_file: str = 'config/genome_browser.yml',
) -> str:
    with open(config_file) as config_handle:
        config = yaml.safe_load(config_handle)
    gene_args = {
        'i': gene_id,
        'c': chromosome,
        'l': left_position,
        'r': right_position,
    }
    return '{base_url}?{gene_args}&{other_args}'.format(
        base_url=config['base_url'],
        gene_args=config_dict_to_get_string(gene_args),
        other_args=config_dict_to_get_string(config['get_args']),
    )


def get_url_contents(url: str, timeout_seconds: int = 60) -> str:
    response = requests.get(url, timeout=timeout_seconds)
    if response.status_code == 200:  # noqa: WPS432
        return response.content.decode('UTF-8')
    raise IOError('Error: {code}'.format(code=response.status_code))


def parse_fasta_string(url_contents: str) -> str:
    soup = BeautifulSoup(url_contents, 'html.parser')
    pre_tag = soup.find('pre')
    if not pre_tag:
        raise ValueError('<pre> tag not found.')
    return pre_tag.text.strip()


def fasta_records_from_url(
    url: str,
    timeout_seconds: int = 60,
) -> list[SeqRecord]:
    url_contents = get_url_contents(url, timeout_seconds)
    fasta_string = parse_fasta_string(url_contents)
    fasta_records: list[SeqRecord] = []
    with io.StringIO(fasta_string) as fasta:
        for fasta_record in SeqIO.parse(fasta, 'fasta'):
            fasta_records.append(fasta_record)
    return fasta_records
