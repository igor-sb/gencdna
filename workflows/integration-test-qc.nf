nextflow.enable.dsl = 2

params.test_file = "$projectDir/../tests/expected_error_filtering/fixtures/reads.fastq"
params.exons = "$projectDir/../tests/exon_alignment/fixtures/exons_14nt.fasta"
params.blast_config = "$projectDir/../tests/exon_alignment/fixtures/blast_config.yml"
// params.test_filtered_file = "$projectDir/../tests/fixtures/"
params.forward_primer = "ATGG"
params.reverse_primer = "GATT"
params.out_dir = ""


workflow {
    filtered_ch = filter_reads_with_low_expected_errors(params.test_file)
    flagged_ch = flag_reads_with_low_quality_repeated_bases(filtered_ch)
    pcr_filtered_ch = filter_reads_without_pcr_primers(
        flagged_ch,
        params.forward_primer,
        params.reverse_primer
    )
    fasta_export_ch = export_unique_reads_to_fasta(pcr_filtered_ch)
    exon_align_ch = align_exons(
        fasta_export_ch,
        params.exons,
        params.blast_config
    )
}


process filter_reads_with_low_expected_errors {
    input:
    path raw_reads

    output:
    path "filtered_reads.fastq.gz"

    script:
    """
    poetry run python \
        $projectDir/../gencdna/file_api/expected_error_filter.py \
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
    poetry run python \
        $projectDir/../gencdna/file_api/repeat_base_flagging.py \
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
    poetry run cutadapt \
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
    path "pcr_filtered_reads.fastq.gz"

    output:
    path "unique_reads.fasta"

    script:
    """
    poetry run python $projectDir/../gencdna/file_api/fastq_to_fasta.py \
        pcr_filtered_reads.fastq.gz \
        unique_reads.fasta
    """
}


process align_exons {
    input:
    path "unique_reads.fasta"
    path exons_fasta
    path blast_config

    output:
    path "blast_output.tsv"

    publishDir "$params.out_dir"

    script:
    """
    poetry run python $projectDir/../gencdna/file_api/exon_alignment.py \
        $exons_fasta \
        unique_reads.fasta \
        blast_output.tsv \
        -b $blast_config 
    """
}