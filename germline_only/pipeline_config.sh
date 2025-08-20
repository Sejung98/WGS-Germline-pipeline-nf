#!/bin/bash

# =============================================================================
# Hereditary Cancer Germline Pipeline Configuration
# =============================================================================

# 기본 디렉토리 설정
export BASE_DIR="/home/ricky8419/09_Hereditary_cancer"
export BAM_DIR="${BASE_DIR}/00_rawdata/bam_markdup"
export OUTPUT_DIR="${BASE_DIR}/Results_germlineMode"
export TOOLS_DIR="/home/ricky8419/HMF_pipeline/tool"
export RESOURCE_DIR="/home/ricky8419/HMF_pipeline/resource/38"

# 참조 유전체 설정
export REF_GENOME_VERSION="38"
export REF_GENOME="${RESOURCE_DIR}/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# 공통 리소스 파일
export ENSEMBL_DIR="${RESOURCE_DIR}/common/ensembl_data"
export DRIVER_GENE_PANEL="${RESOURCE_DIR}/common/DriverGenePanel.38.tsv"

# SAGE & PAVE 관련 파일
export SAGE_JAR="${TOOLS_DIR}/sage.jar"
export PAVE_JAR="${TOOLS_DIR}/pave.jar"
export HOTSPOTS="${RESOURCE_DIR}/variants/KnownHotspots.germline.38.vcf.gz"
export PANEL_BED="${RESOURCE_DIR}/variants/ActionableCodingPanel.38.bed.gz"
export COVERAGE_BED="${RESOURCE_DIR}/variants/CoverageCodingPanel.38.bed.gz"
export HIGH_CONFIDENCE_BED="${RESOURCE_DIR}/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz"
export MAPPABILITY_BED="${RESOURCE_DIR}/variants/mappability_150.38.bed.gz"
export CLINVAR_VCF="${RESOURCE_DIR}/variants/clinvar.38.vcf.gz"
export BLACKLIST_BED="${RESOURCE_DIR}/variants/KnownBlacklist.germline.38.bed"
export BLACKLIST_VCF="${RESOURCE_DIR}/variants/KnownBlacklist.germline.38.vcf.gz"

# AMBER & COBALT 관련 파일
export AMBER_JAR="${TOOLS_DIR}/amber.jar"
export COBALT_JAR="${TOOLS_DIR}/cobalt.jar"
export LOCI_FILE="${RESOURCE_DIR}/copy_number/GermlineHetPon.38.vcf.gz"
export GC_PROFILE="${RESOURCE_DIR}/copy_number/GC_profile.1000bp.38.cnp"

# GRIDSS 관련 파일
export SV_PREP_JAR="${TOOLS_DIR}/sv-prep.jar"
export GRIDSS_SCRIPT="${TOOLS_DIR}/gridss.run.sh"
export GRIDSS_JAR="${TOOLS_DIR}/gridss.jar"
export SV_BLACKLIST_BED="${RESOURCE_DIR}/sv/sv_prep_blacklist.38.bed"
export KNOWN_FUSION_BED="${RESOURCE_DIR}/sv/known_fusions.38.bedpe"
export GRIDSS_BLACKLIST_BED="${RESOURCE_DIR}/sv/gridss_blacklist.38.bed"
export GRIDSS_CONFIG="${RESOURCE_DIR}/sv/gridss.properties"

# PURPLE 관련 파일
export PURPLE_JAR="${TOOLS_DIR}/purple.jar"
export GERMLINE_HOTSPOTS="${RESOURCE_DIR}/variants/KnownHotspots.germline.38.vcf.gz"
export GERMLINE_DEL_FREQ_FILE="${RESOURCE_DIR}/copy_number/cohort_germline_del_freq.38.csv"

# GRIPSS 관련 파일
export GRIPSS_JAR="${TOOLS_DIR}/gripss.jar"
export SV_HOTSPOT="${RESOURCE_DIR}/sv/known_fusions.38.bedpe"
export SV_PON_FILE="${RESOURCE_DIR}/sv/sv_pon.38.bedpe.gz"
export SGL_PON_FILE="${RESOURCE_DIR}/sv/sgl_pon.38.bed.gz"
export REPEAT_MASK_FILE="${RESOURCE_DIR}/sv/repeat_mask_data.38.fa.gz"

# LINX 관련 파일
export LINX_JAR="${TOOLS_DIR}/linx.jar"

# 시스템 리소스 설정
export THREADS=64
export MAX_MEMORY=512
export GRIDSS_THREADS=8

# 로그 설정
export LOG_DIR="${OUTPUT_DIR}/logs"
export TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# 출력 디렉토리 생성
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

echo "Configuration loaded successfully!"
echo "Base directory: ${BASE_DIR}"
echo "BAM directory: ${BAM_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo "Memory: ${MAX_MEMORY}G"
