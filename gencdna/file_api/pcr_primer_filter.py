"""CLI for PCR primer filtering from a FASTA/FASTQ(.GZ) with cutadapt."""

import logging

import fire

from gencdna.bash import BinaryExecWithYamlArgs

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


def filter_reads_without_pcr_primers(
    input_file: str,
    output_file: str,
    forward_primer: str,
    reverse_primer: str,
    config_yml: str = 'config/cutadapt.yml',
) -> str:
    cutadapt = BinaryExecWithYamlArgs('cutadapt', config_yml)
    cutadapt.config['arguments'] = {
        '-g': '{forward_primer}...{reverse_primer}'.format(
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
        ),
        **cutadapt.config['arguments'],
    }
    if output_file:
        cutadapt.config['arguments']['--output'] = output_file
    cutadapt_output = cutadapt.run(input_file)
    if cutadapt_output.stderr != b'':
        raise RuntimeError(cutadapt_output.stderr.decode('UTF-8'))
    return cutadapt_output.stdout.decode('UTF-8')


if __name__ == '__main__':
    fire.Fire(filter_reads_without_pcr_primers)
