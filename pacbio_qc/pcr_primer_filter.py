"""PCR primer filtering from a FASTA/FASTQ file using cutadapt."""

from pacbio_qc.bash import BinaryExecWithYamlArgs


def find_reads_with_pcr_primers(
    input_fastx_file: str,
    cutadapt_args: dict[str, str],
) -> str:
    cutadapt = BinaryExecWithYamlArgs('cutadapt')
    cutadapt.config['arguments'] = cutadapt_args
    cutadapt_output = cutadapt.run(input_fastx_file)
    if cutadapt_output.stderr != b'':
        raise RuntimeError(cutadapt_output.stderr.decode('UTF-8'))
    return cutadapt_output.stdout.decode('UTF-8')
