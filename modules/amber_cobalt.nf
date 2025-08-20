process AMBER_COBALT {
    tag "$sample_id"
    
    cpus params.threads
    memory params.memory
    
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.amber.baf.tsv.gz"), path("${sample_id}.cobalt.ratio.tsv.gz")
    
    when:
    bam.exists()
    
    script:
    """
    echo "=== AMBER START: \$(date) ===" 
    
    java -Xmx${task.memory.toGiga()}G \\
        -jar ${params.amber_jar} \\
        -reference ${sample_id} \\
        -reference_bam ${bam} \\
        -ref_genome_version ${params.ref_genome_version} \\
        -output_dir . \\
        -threads ${task.cpus} \\
        -loci ${params.loci_file}
    
    echo "=== AMBER END: \$(date) ==="
    echo "=== COBALT START: \$(date) ==="
    
    java -Xmx${task.memory.toGiga()}G \\
        -jar ${params.cobalt_jar} \\
        -reference ${sample_id} \\
        -reference_bam ${bam} \\
        -output_dir . \\
        -threads ${task.cpus} \\
        -gc_profile ${params.gc_profile}
    
    echo "=== COBALT END: \$(date) ==="
    """
}
