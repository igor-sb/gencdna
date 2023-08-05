from pacbio_qc.bash import BinaryExecWithYamlArgs


def filter_reads_without_pcr_primers(
    input_fasta_file: str,
    output_fasta_file: str
) -> None:
    cutadapt = BinaryExecWithYamlArgs('cutadapt', )
