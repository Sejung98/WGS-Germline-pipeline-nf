process GRIDSS {
    tag "$sample_id"
    
    cpus params.threads
    memory params.memory
    
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.gridss.raw.vcf.gz"), path("${sample_id}_gridss.vcf.gz")
    
    when:
    bam.exists()
    
    script:
    """
    echo "=== SvPrep START: \$(date) ===" 
    
    java -Xmx${task.memory.toGiga()}G -jar ${params.sv_prep_jar} \\
        -sample ${sample_id} \\
        -bam_file ${bam} \\
        -ref_genome ${params.ref_genome} \\
        -ref_genome_version ${params.ref_genome_version} \\
        -blacklist_bed ${params.sv_blacklist_bed} \\
        -known_fusion_bed ${params.known_fusion_bed} \\
        -write_types "JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST" \\
        -output_dir . \\
        -threads ${task.cpus}
    
    echo "=== SvPrep END: \$(date) ==="
    echo "=== BAM Sort START: \$(date) ==="
    
    # Sort and index the sv_prep BAM
    samtools sort -@ ${task.cpus} -m 2G -T tmp_${sample_id} -O bam ${sample_id}.sv_prep.bam -o ${sample_id}.sv_prep.sorted.bam
    samtools index -@ ${task.cpus} ${sample_id}.sv_prep.sorted.bam
    rm -f ${sample_id}.sv_prep.bam
    
    echo "=== BAM Sort END: \$(date) ==="
    echo "=== GRIDSS START: \$(date) ==="
    
    # Run GRIDSS
    ${params.gridss_script} --jar ${params.gridss_jar} \\
        --steps all \\
        --workingdir . \\
        --reference ${params.ref_genome} \\
        --blacklist ${params.gridss_blacklist_bed} \\
        --configuration ${params.gridss_config} \\
        --labels ${sample_id} \\
        --bams ${bam} \\
        --filtered_bams ${sample_id}.sv_prep.sorted.bam \\
        --output ${sample_id}.gridss.raw.vcf.gz \\
        --jvmheap ${task.memory.toGiga()}G \\
        --threads 8
    
    echo "=== GRIDSS END: \$(date) ==="
    echo "=== DepthAnnotator START: \$(date) ==="
    
    # Run SvPrep Reference Depth Annotator
    java -Xmx${task.memory.toGiga()}G -cp ${params.sv_prep_jar} com.hartwig.hmftools.svprep.depth.DepthAnnotator \\
        -input_vcf ${sample_id}.gridss.raw.vcf.gz \\
        -output_vcf ${sample_id}_gridss.vcf.gz \\
        -samples ${sample_id} \\
        -bam_files ${bam} \\
        -ref_genome ${params.ref_genome} \\
        -ref_genome_version ${params.ref_genome_version} \\
        -threads ${task.cpus}
    
    echo "=== DepthAnnotator END: \$(date) ==="
    
    # Cleanup
    rm -rf tmp_${sample_id}*
    """
}
