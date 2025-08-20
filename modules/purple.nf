process PURPLE {
    tag "$sample_id"
    
    cpus params.threads
    memory params.memory
    
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(pave_vcf), path(gripss_vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.purple.sv.germline.vcf.gz")
    
    when:
    pave_vcf.exists() && gripss_vcf.exists()
    
    script:
    """
    echo "=== PURPLE START: \$(date) ===" 
    
    java -Xmx${task.memory.toGiga()}G \\
        -jar ${params.purple_jar} \\
        -reference ${sample_id} \\
        -amber ${params.output_dir} \\
        -cobalt ${params.output_dir} \\
        -gc_profile ${params.gc_profile} \\
        -ref_genome_version ${params.ref_genome_version} \\
        -ref_genome ${params.ref_genome} \\
        -ensembl_data_dir ${params.ensembl_dir} \\
        -germline_hotspots ${params.germline_hotspots} \\
        -germline_del_freq_file ${params.germline_del_freq_file} \\
        -germline_vcf ${pave_vcf} \\
        -germline_sv_vcf ${gripss_vcf} \\
        -driver_gene_panel ${params.driver_gene_panel} \\
        -output_dir . \\
        -threads ${task.cpus}
    
    echo "=== PURPLE END: \$(date) ==="
    """
}
