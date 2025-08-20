process GRIPSS {
    tag "$sample_id"
    
    cpus 4
    memory params.memory
    
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(gridss_raw_vcf), path(gridss_vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.gripss.filtered.germline.vcf.gz")
    
    when:
    gridss_raw_vcf.exists()
    
    script:
    """
    echo "=== GRIPSS START: \$(date) ===" 
    
    java -Xmx${task.memory.toGiga()}G -jar ${params.gripss_jar} \\
        -sample ${sample_id} \\
        -germline \\
        -ref_genome ${params.ref_genome} \\
        -ref_genome_version ${params.ref_genome_version} \\
        -known_hotspot_file ${params.sv_hotspot} \\
        -pon_sv_file ${params.sv_pon_file} \\
        -pon_sgl_file ${params.sgl_pon_file} \\
        -repeat_mask_file ${params.repeat_mask_file} \\
        -vcf ${gridss_raw_vcf} \\
        -output_dir . \\
        -output_id germline
    
    echo "=== GRIPSS END: \$(date) ==="
    """
}
