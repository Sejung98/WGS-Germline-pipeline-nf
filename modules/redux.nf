#!/usr/bin/env nextflow

process REDUX {
    tag "${sample_id}"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.redux.bam"), path("${sample_id}.jitter_params.tsv"), path("${sample_id}.ms_table.tsv.gz"), emit: redux_output
    path "*.log", optional: true
    
    script:
    def redux_jar = "${params.tools_dir}/redux.jar"
    def unmap_regions = "${params.resource_dir}/common/unmap_regions.38.tsv"
    def msi_sites = "${params.resource_dir}/common/msi_jitter_sites.38.tsv.gz"
    def samtools = "${params.tools_dir}/samtools"
    
    """
    #!/bin/bash
    set -e
    
    echo "=== REDUX START: \$(date) ==="
    echo "Sample: ${sample_id}"
    echo "Input BAM: ${bam_file}"
    
    # REDUX 실행
    java -Xmx${params.memory} -jar ${redux_jar} \\
        -sample "${sample_id}" \\
        -input_bam "${bam_file}" \\
        -output_bam "${sample_id}.redux.bam" \\
        -ref_genome "${params.ref_genome}" \\
        -ref_genome_version "V${params.ref_genome_version}" \\
        -unmap_regions "${unmap_regions}" \\
        -ref_genome_msi_file "${msi_sites}" \\
        -form_consensus \\
        -write_stats \\
        -bamtool "${samtools}" \\
        -output_dir "." \\
        -threads ${params.threads} \\
        -log_level INFO
    
    # 출력 파일 검증
    if [[ ! -f "${sample_id}.redux.bam" ]]; then
        echo "ERROR: REDUX BAM file not created"
        exit 1
    fi
    
    if [[ ! -f "${sample_id}.jitter_params.tsv" ]]; then
        echo "ERROR: Jitter parameters file not created"
        exit 1
    fi
    
    if [[ ! -f "${sample_id}.ms_table.tsv.gz" ]]; then
        echo "ERROR: MS table file not created"
        exit 1
    fi
    
    # BAM 인덱스 생성
    ${samtools} index "${sample_id}.redux.bam"
    
    echo "=== REDUX END: \$(date) ==="
    echo "Output files:"
    echo "  - ${sample_id}.redux.bam"
    echo "  - ${sample_id}.jitter_params.tsv"
    echo "  - ${sample_id}.ms_table.tsv.gz"
    """
}
