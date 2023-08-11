nextflow.enable.dsl = 2

params.test_file = "$projectDir/../tests/fixtures/example_reads_with_exons.fastq.gz"
params.test_filtered_file = "$projectDir/../tests/fixtures/"
params.forward_primer = "ATGG"
params.reverse_primer = "GATT"


process filter_reads_with_low_expected_errors {
    input:
    path raw_reads

    output:
    path "filtered_reads.fastq.gz"

    script:
    """
    python \
        $projectDir/../pacbio_qc/file_api/expected_error_filter.py \
        $raw_reads \
        filtered_reads.fastq.gz \
        0.01
    """
}

process flag_reads_with_low_quality_repeated_bases {
    input:
    path "filtered_reads.fastq.gz"

    output:
    path "flagged_reads.fastq.gz"

    script:
    """
    python \
        $projectDir/../pacbio_qc/file_api/low_quality_repeated_base_flagger.py \
        filtered_reads.fastq.gz \
        flagged_reads.fastq.gz
    """
}


process filter_reads_without_pcr_primers {
    input:
    path "flagged_reads.fastq.gz"
    val forward_primer
    val reverse_primer

    output:
    path "pcr_filtered_reads.fastq.gz"

    script:
    """
    cutadapt \
        -g $forward_primer...$reverse_primer \
        --trimmed-only \
        --minimum-length 1 \
        --revcomp \
        --quiet \
        --output pcr_filtered_reads.fastq.gz \
        flagged_reads.fastq.gz
    """
}


process export_unique_reads_to_fasta {
    input:
    "pcr_filtered_reads.fastq.gz"

    output:
    "unique_reads.fasta"

    script:
    """
    """
}


workflow {
    filtered_ch = filter_reads_with_low_expected_errors(params.test_file)
    flagged_ch = flag_reads_with_low_quality_repeated_bases(filtered_ch)
    pcr_filtered_ch = filter_reads_without_pcr_primers(
        flagged_ch,
        params.forward_primer,
        params.reverse_primer
    )
}
