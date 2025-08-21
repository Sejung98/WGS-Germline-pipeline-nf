#!/usr/bin/env nextflow

/*
 * Hereditary Cancer Germline Analysis Pipeline
 * Nextflow implementation
 */

nextflow.enable.dsl=2

// Pipeline parameters
params.bam_dir = "/home/ricky8419/09_Hereditary_cancer/00_rawdata/bam_markdup"
params.output_dir = "/home/ricky8419/09_Hereditary_cancer/Results_germlineMode_nextflow"
params.tools_dir = "/home/ricky8419/HMF_pipeline/tool"
params.resource_dir = "/home/ricky8419/HMF_pipeline/resource/38"
params.ref_genome_version = "38"
params.threads = 8
params.memory = "64.GB"

// Tool paths
params.sage_jar = "${params.tools_dir}/sage.jar"
params.pave_jar = "${params.tools_dir}/pave.jar"
params.amber_jar = "${params.tools_dir}/amber.jar"
params.cobalt_jar = "${params.tools_dir}/cobalt.jar"
params.sv_prep_jar = "${params.tools_dir}/sv-prep.jar"
params.gridss_jar = "${params.tools_dir}/gridss.jar"
params.gridss_script = "${params.tools_dir}/gridss.run.sh"
params.purple_jar = "${params.tools_dir}/purple.jar"
params.gripss_jar = "${params.tools_dir}/gripss.jar"
params.linx_jar = "${params.tools_dir}/linx.jar"
params.redux_jar = "${params.tools_dir}/redux.jar"
params.samtools = "${params.tools_dir}/samtools"

// Reference files
params.ref_genome = "${params.resource_dir}/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
params.ensembl_dir = "${params.resource_dir}/common/ensembl_data"
params.driver_gene_panel = "${params.resource_dir}/common/DriverGenePanel.38.tsv"
params.hotspots = "${params.resource_dir}/variants/KnownHotspots.germline.38.vcf.gz"
params.high_confidence_bed = "${params.resource_dir}/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz"
params.mappability_bed = "${params.resource_dir}/variants/mappability_150.38.bed.gz"
params.clinvar_vcf = "${params.resource_dir}/variants/clinvar.38.vcf.gz"
params.blacklist_bed = "${params.resource_dir}/variants/KnownBlacklist.germline.38.bed"
params.blacklist_vcf = "${params.resource_dir}/variants/KnownBlacklist.germline.38.vcf.gz"

// AMBER & COBALT
params.loci_file = "${params.resource_dir}/copy_number/GermlineHetPon.38.vcf.gz"
params.gc_profile = "${params.resource_dir}/copy_number/GC_profile.1000bp.38.cnp"

// GRIDSS
params.sv_blacklist_bed = "${params.resource_dir}/sv/sv_prep_blacklist.38.bed"
params.known_fusion_bed = "${params.resource_dir}/sv/known_fusions.38.bedpe"
params.gridss_blacklist_bed = "${params.resource_dir}/sv/gridss_blacklist.38.bed"
params.gridss_config = "${params.resource_dir}/sv/gridss.properties"

// PURPLE
params.germline_hotspots = "${params.resource_dir}/variants/KnownHotspots.germline.38.vcf.gz"
params.germline_del_freq_file = "${params.resource_dir}/copy_number/cohort_germline_del_freq.38.csv"

// GRIPSS
params.sv_hotspot = "${params.resource_dir}/sv/known_fusions.38.bedpe"
params.sv_pon_file = "${params.resource_dir}/sv/sv_pon.38.bedpe.gz"
params.sgl_pon_file = "${params.resource_dir}/sv/sgl_pon.38.bed.gz"
params.repeat_mask_file = "${params.resource_dir}/sv/repeat_mask_data.38.fa.gz"

log.info """\
==============================================
H E R E D I T A R Y   C A N C E R   P I P E L I N E
==============================================
BAM directory     : ${params.bam_dir}
Output directory  : ${params.output_dir}
Tools directory   : ${params.tools_dir}
Resource directory: ${params.resource_dir}
Threads per job   : ${params.threads}
Memory per job    : ${params.memory}
==============================================
"""

// Import modules
include { REDUX } from './modules/redux'
include { SAGE_PAVE } from './modules/sage_pave'
include { AMBER_COBALT } from './modules/amber_cobalt'
include { GRIDSS } from './modules/gridss'
include { GRIPSS } from './modules/gripss'
include { PURPLE } from './modules/purple'
include { LINX } from './modules/linx'

workflow {
    // Create input channel from BAM files
    bam_ch = Channel
        .fromPath("${params.bam_dir}/*_markdup.bam")
        .map { bam ->
            def sample_id = bam.baseName.replaceAll('_markdup$', '')
            tuple(sample_id, bam)
        }

    // Stage 1: REDUX (BAM 정제 및 jitter 모델링)
    redux_ch = REDUX(bam_ch)

    // Stage 2: SAGE_PAVE (REDUX BAM 사용, jitter 파라미터 적용)
    sage_pave_ch = SAGE_PAVE(redux_ch)

    // Stage 3: AMBER_COBALT (원본 BAM 사용)
    amber_cobalt_ch = AMBER_COBALT(bam_ch)

    // Stage 4: GRIDSS (원본 BAM 사용)
    gridss_ch = GRIDSS(bam_ch)

    // Stage 5: GRIPSS (depends on GRIDSS)
    gripss_ch = GRIPSS(gridss_ch)

    // Combine results for PURPLE
    purple_input_ch = sage_pave_ch
        .join(amber_cobalt_ch)
        .join(gripss_ch)
        .map { sample_id, sage_pave, amber_cobalt, gripss ->
            tuple(sample_id, sage_pave[1], gripss[1])  // pave_vcf, gripss_vcf
        }

    // Stage 6: PURPLE (depends on SAGE_PAVE, AMBER_COBALT, GRIPSS)
    purple_ch = PURPLE(purple_input_ch)

    // Stage 7: LINX (depends on PURPLE)
    linx_ch = LINX(purple_ch)
}

workflow.onComplete {
    log.info """\
==============================================
Pipeline completed at: ${workflow.complete}
Execution status: ${workflow.success ? 'OK' : 'failed'}
Results directory: ${params.output_dir}
==============================================
    """
}
