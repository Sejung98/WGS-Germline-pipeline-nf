process SAGE_PAVE {
    tag "$sample_id"
    
    cpus params.threads
    memory params.memory
    
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.sage.germline.vcf.gz"), path("${sample_id}.pave.germline.vcf.gz")
    
    when:
    bam.exists()
    
    script:
    """
    echo "=== SAGE START: \$(date) ===" 
    
    java -Xmx${task.memory.toGiga()}G -jar ${params.sage_jar} \\
        -tumor ${sample_id} \\
        -tumor_bam ${bam} \\
        -hotspots ${params.hotspots} \\
        -high_confidence_bed ${params.high_confidence_bed} \\
        -ref_genome ${params.ref_genome} \\
        -ref_genome_version ${params.ref_genome_version} \\
        -ensembl_data_dir ${params.ensembl_dir} \\
        -ref_sample_count 0 \\
        -germline \\
        -threads ${task.cpus} \\
        -output_vcf ${sample_id}.sage.germline.vcf.gz
    
    echo "=== SAGE END: \$(date) ==="
    echo "=== PAVE START: \$(date) ===" 
    
    java -Xmx${task.memory.toGiga()}G -jar ${params.pave_jar} \\
        -sample ${sample_id} \\
        -vcf_file ${sample_id}.sage.germline.vcf.gz \\
        -ensembl_data_dir ${params.ensembl_dir} \\
        -driver_gene_panel ${params.driver_gene_panel} \\
        -ref_genome ${params.ref_genome} \\
        -ref_genome_version ${params.ref_genome_version} \\
        -mappability_bed ${params.mappability_bed} \\
        -clinvar_vcf ${params.clinvar_vcf} \\
        -blacklist_bed ${params.blacklist_bed} \\
        -blacklist_vcf ${params.blacklist_vcf} \\
        -read_pass_only \\
        -output_dir . \\
        -output_vcf_file ${sample_id}.pave.germline.vcf.gz
    
    echo "=== PAVE END: \$(date) ==="
    """
}
