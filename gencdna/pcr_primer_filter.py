"""PCR primer filtering from a FASTA/FASTQ file using cutadapt."""

from gencdna.external_calls import BinaryExecWithYamlArgs


def find_reads_with_pcr_primers(
    fastx_file: str,
    cutadapt_args: dict[str, str],
) -> str:
    """Load FASTA/FASTQ file and keep only reads with PCR primers.

    Standard output is returned as a string. To suppress summary statistics,
    add --quiet in cutadapt_args. To output sequences in FASTA/FASTQ, add
    --output filename.fasta to cutadapt_args.

    Raises:
        RuntimeError: If cutadapt throws a standard error.
    """
    cutadapt = BinaryExecWithYamlArgs('cutadapt')
    cutadapt.config['arguments'] = cutadapt_args
    cutadapt_output = cutadapt.run(fastx_file)
    if cutadapt_output.stderr != b'':
        raise RuntimeError(cutadapt_output.stderr.decode('UTF-8'))
    return cutadapt_output.stdout.decode('UTF-8')
