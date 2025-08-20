process LINX {
    tag "$sample_id"
    
    cpus 4
    memory params.memory
    
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(purple_sv_vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.linx.germline.disruption.tsv"), path("${sample_id}.linx.germline.breakend.tsv")
    
    when:
    purple_sv_vcf.exists()
    
    script:
    """
    echo "=== LINX START: \$(date) ===" 
    
    java -Xmx${task.memory.toGiga()}G \\
        -jar ${params.linx_jar} \\
        -sample ${sample_id} \\
        -germline \\
        -ref_genome_version ${params.ref_genome_version} \\
        -sv_vcf ${purple_sv_vcf} \\
        -output_dir . \\
        -ensembl_data_dir ${params.ensembl_dir} \\
        -driver_gene_panel ${params.driver_gene_panel}
    
    echo "=== LINX END: \$(date) ==="
    """
}
