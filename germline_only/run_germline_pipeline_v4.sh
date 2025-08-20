#!/bin/bash

# =============================================================================
# Hereditary Cancer Germline Pipeline Main Script v4
# Each sample runs its own complete pipeline independently
# =============================================================================

set -e  # 오류 시 스크립트 중단

# 스크립트 디렉토리 경로
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 설정 파일 로드
source "${SCRIPT_DIR}/pipeline_config.sh"

# 전역 로그 파일
MAIN_LOG="${LOG_DIR}/pipeline_main_${TIMESTAMP}.log"

# 로그 함수들
log() {
    local message="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "${message}" | tee -a "${MAIN_LOG}"
}

log_info() {
    log "INFO: $1"
}

log_warn() {
    log "WARN: $1"
}

log_error() {
    log "ERROR: $1"
}

error_exit() {
    log_error "$1"
    exit 1
}

# 사용법 출력
show_usage() {
    cat << EOF
Usage: $0 [OPTIONS] [SAMPLE_ID1] [SAMPLE_ID2] ...

이 스크립트는 각 샘플별로 독립적인 전체 파이프라인을 실행합니다.
각 샘플은 1단계→2단계→3단계 순서로 완전히 처리됩니다.

OPTIONS:
  -h, --help              Show this help message
  -l, --list              List available samples
  -a, --all               Process all samples (default)
  --parallel JOBS         Number of parallel sample processes (default: 4)
  --dry-run               Show what would be executed without running
  --max-memory MEM        Memory per sample process in GB (default: from config)
  --threads NUM           Threads per sample process (default: from config)

EXAMPLES:
  $0                      # Process all samples, each with full pipeline
  $0 sample1 sample2      # Process specific samples only
  $0 --parallel 8         # Run 8 samples in parallel
  $0 --list               # List available samples
  $0 --dry-run            # Show execution plan
  $0 --max-memory 128     # Use 128GB per sample
  $0 --threads 16         # Use 16 threads per sample

PIPELINE STAGES (each sample runs all stages):
  Stage 1: SAGE_PAVE, AMBER_COBALT, GRIDSS (parallel within sample)
  Stage 2: GRIPSS, PURPLE (sequential, depends on Stage 1)
  Stage 3: LINX (depends on Stage 2)
EOF
}

# BAM 파일 존재 확인
check_bam_files() {
    log_info "Checking BAM files in ${BAM_DIR}..."
    if [[ ! -d "${BAM_DIR}" ]]; then
        error_exit "BAM directory does not exist: ${BAM_DIR}"
    fi
    
    bam_count=$(find "${BAM_DIR}" -name "*_markdup.bam" | wc -l)
    if [[ ${bam_count} -eq 0 ]]; then
        error_exit "No _markdup.bam files found in ${BAM_DIR}"
    fi
    
    log_info "Found ${bam_count} BAM files"
}

# 사용 가능한 샘플 ID 목록 생성
get_available_samples() {
    available_samples=()
    for bam_file in "${BAM_DIR}"/*_markdup.bam; do
        if [[ -f "${bam_file}" ]]; then
            sample_id=$(basename "${bam_file}" "_markdup.bam")
            available_samples+=("${sample_id}")
        fi
    done
    
    if [[ ${#available_samples[@]} -eq 0 ]]; then
        error_exit "No valid BAM files found"
    fi
}

# 샘플 목록 출력
list_samples() {
    get_available_samples
    log_info "Available samples (${#available_samples[@]}):"
    for sample_id in "${available_samples[@]}"; do
        echo "  - ${sample_id}"
    done
}

# 샘플 유효성 검사
validate_samples() {
    local requested_samples=("$@")
    local invalid_samples=()
    
    get_available_samples
    
    for requested in "${requested_samples[@]}"; do
        local found=false
        for available in "${available_samples[@]}"; do
            if [[ "${requested}" == "${available}" ]]; then
                found=true
                break
            fi
        done
        if [[ "${found}" == "false" ]]; then
            invalid_samples+=("${requested}")
        fi
    done
    
    if [[ ${#invalid_samples[@]} -gt 0 ]]; then
        log_error "Invalid sample IDs: ${invalid_samples[*]}"
        log_info "Available samples: ${available_samples[*]}"
        exit 1
    fi
}

# SAGE & PAVE 실행 함수
run_sage_pave() {
    local sample_id=$1
    local sample_memory=$2
    local sample_threads=$3
    local log_file="${LOG_DIR}/sage_pave_${sample_id}_${TIMESTAMP}.log"
    
    log_info "[${sample_id}] Starting SAGE & PAVE"
    
    local reference_bam="${BAM_DIR}/${sample_id}_markdup.bam"
    local sage_output_vcf="${OUTPUT_DIR}/${sample_id}.sage.germline.vcf.gz"
    local pave_output_vcf="${OUTPUT_DIR}/${sample_id}.pave.germline.vcf.gz"
    
    # 입력 파일 확인
    if [[ ! -f "${reference_bam}" ]]; then
        error_exit "[${sample_id}] BAM file not found: ${reference_bam}"
    fi
    
    # SAGE 실행
    if [[ -f "${sage_output_vcf}" ]]; then
        log_info "[${sample_id}] Skipping SAGE, output exists"
    else
        log_info "[${sample_id}] Running SAGE"
        echo "=== SAGE START: $(date) ===" >> "${log_file}"
        
        java -Xmx${sample_memory}G -jar "${SAGE_JAR}" \
            -tumor "${sample_id}" \
            -tumor_bam "${reference_bam}" \
            -hotspots "${HOTSPOTS}" \
            -high_confidence_bed "${HIGH_CONFIDENCE_BED}" \
            -ref_genome "${REF_GENOME}" \
            -ref_genome_version "${REF_GENOME_VERSION}" \
            -ensembl_data_dir "${ENSEMBL_DIR}" \
            -ref_sample_count 0 \
            -germline \
            -threads "${sample_threads}" \
            -output_vcf "${sage_output_vcf}" >> "${log_file}" 2>&1
        
        local sage_exit_code=$?
        echo "=== SAGE END: $(date), Exit Code: ${sage_exit_code} ===" >> "${log_file}"
        
        if [[ ${sage_exit_code} -eq 0 ]]; then
            log_info "[${sample_id}] SAGE completed successfully"
        else
            error_exit "[${sample_id}] SAGE failed. Check log: ${log_file}"
        fi
    fi
    
    # PAVE 실행
    if [[ -f "${pave_output_vcf}" ]]; then
        log_info "[${sample_id}] Skipping PAVE, output exists"
    else
        log_info "[${sample_id}] Running PAVE"
        echo "=== PAVE START: $(date) ===" >> "${log_file}"
        
        java -Xmx${sample_memory}G -jar "${PAVE_JAR}" \
            -sample "${sample_id}" \
            -vcf_file "${sage_output_vcf}" \
            -ensembl_data_dir "${ENSEMBL_DIR}" \
            -driver_gene_panel "${DRIVER_GENE_PANEL}" \
            -ref_genome "${REF_GENOME}" \
            -ref_genome_version "${REF_GENOME_VERSION}" \
            -mappability_bed "${MAPPABILITY_BED}" \
            -clinvar_vcf "${CLINVAR_VCF}" \
            -blacklist_bed "${BLACKLIST_BED}" \
            -blacklist_vcf "${BLACKLIST_VCF}" \
            -read_pass_only \
            -output_dir "${OUTPUT_DIR}" \
            -output_vcf_file "${pave_output_vcf}" >> "${log_file}" 2>&1
        
        local pave_exit_code=$?
        echo "=== PAVE END: $(date), Exit Code: ${pave_exit_code} ===" >> "${log_file}"
        
        if [[ ${pave_exit_code} -eq 0 ]]; then
            log_info "[${sample_id}] PAVE completed successfully"
        else
            error_exit "[${sample_id}] PAVE failed. Check log: ${log_file}"
        fi
    fi
}

# AMBER & COBALT 실행 함수
run_amber_cobalt() {
    local sample_id=$1
    local sample_memory=$2
    local sample_threads=$3
    local log_file="${LOG_DIR}/amber_cobalt_${sample_id}_${TIMESTAMP}.log"
    
    log_info "[${sample_id}] Starting AMBER & COBALT"
    
    local reference_bam="${BAM_DIR}/${sample_id}_markdup.bam"
    
    # 입력 파일 확인
    if [[ ! -f "${reference_bam}" ]]; then
        error_exit "[${sample_id}] BAM file not found: ${reference_bam}"
    fi
    
    # AMBER 실행
    log_info "[${sample_id}] Running AMBER"
    echo "=== AMBER START: $(date) ===" >> "${log_file}"
    
    java -Xmx${sample_memory}G \
        -jar "${AMBER_JAR}" \
        -reference "${sample_id}" \
        -reference_bam "${reference_bam}" \
        -ref_genome_version "${REF_GENOME_VERSION}" \
        -output_dir "${OUTPUT_DIR}" \
        -threads "${sample_threads}" \
        -loci "${LOCI_FILE}" >> "${log_file}" 2>&1
    
    local amber_exit_code=$?
    echo "=== AMBER END: $(date), Exit Code: ${amber_exit_code} ===" >> "${log_file}"
    
    if [[ ${amber_exit_code} -ne 0 ]]; then
        error_exit "[${sample_id}] AMBER failed. Check log: ${log_file}"
    fi
    
    # COBALT 실행
    log_info "[${sample_id}] Running COBALT"
    echo "=== COBALT START: $(date) ===" >> "${log_file}"
    
    java -Xmx${sample_memory}G \
        -jar "${COBALT_JAR}" \
        -reference "${sample_id}" \
        -reference_bam "${reference_bam}" \
        -output_dir "${OUTPUT_DIR}" \
        -threads "${sample_threads}" \
        -gc_profile "${GC_PROFILE}" >> "${log_file}" 2>&1
    
    local cobalt_exit_code=$?
    echo "=== COBALT END: $(date), Exit Code: ${cobalt_exit_code} ===" >> "${log_file}"
    
    if [[ ${cobalt_exit_code} -ne 0 ]]; then
        error_exit "[${sample_id}] COBALT failed. Check log: ${log_file}"
    fi
    
    log_info "[${sample_id}] AMBER & COBALT completed successfully"
}

# GRIDSS 실행 함수
run_gridss() {
    local sample_id=$1
    local sample_memory=$2
    local sample_threads=$3
    local log_file="${LOG_DIR}/gridss_${sample_id}_${TIMESTAMP}.log"
    
    log_info "[${sample_id}] Starting GRIDSS"
    
    local reference_bam="${BAM_DIR}/${sample_id}_markdup.bam"
    local sv_prep_ref_bam="${OUTPUT_DIR}/${sample_id}.sv_prep.bam"
    local sv_prep_ref_sorted_bam="${OUTPUT_DIR}/${sample_id}.sv_prep.sorted.bam"
    local gridss_raw_vcf="${OUTPUT_DIR}/${sample_id}.gridss.raw.vcf.gz"
    local gridss_vcf="${OUTPUT_DIR}/${sample_id}_gridss.vcf.gz"
    
    # 입력 파일 확인
    if [[ ! -f "${reference_bam}" ]]; then
        error_exit "[${sample_id}] BAM file not found: ${reference_bam}"
    fi
    
    # SvPrep 실행
    local sv_prep_write_types="JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST"
    
    log_info "[${sample_id}] Running SvPrep"
    echo "=== SvPrep START: $(date) ===" >> "${log_file}"
    
    java -Xmx${sample_memory}G -jar "${SV_PREP_JAR}" \
        -sample "${sample_id}" \
        -bam_file "${reference_bam}" \
        -ref_genome "${REF_GENOME}" \
        -ref_genome_version "${REF_GENOME_VERSION}" \
        -blacklist_bed "${SV_BLACKLIST_BED}" \
        -known_fusion_bed "${KNOWN_FUSION_BED}" \
        -write_types "${sv_prep_write_types}" \
        -output_dir "${OUTPUT_DIR}" \
        -threads "${sample_threads}" >> "${log_file}" 2>&1
    
    local svprep_exit_code=$?
    echo "=== SvPrep END: $(date), Exit Code: ${svprep_exit_code} ===" >> "${log_file}"
    
    if [[ ${svprep_exit_code} -ne 0 ]]; then
        error_exit "[${sample_id}] SvPrep failed. Check log: ${log_file}"
    fi
    
    # BAM 파일 정렬 및 인덱싱
    log_info "[${sample_id}] Sorting and indexing BAM"
    echo "=== BAM Sort START: $(date) ===" >> "${log_file}"
    
    samtools sort -@ "${sample_threads}" -m 2G -T "tmp_${sample_id}" -O bam "${sv_prep_ref_bam}" -o "${sv_prep_ref_sorted_bam}" >> "${log_file}" 2>&1
    samtools index -@ "${sample_threads}" "${sv_prep_ref_sorted_bam}" >> "${log_file}" 2>&1
    
    echo "=== BAM Sort END: $(date) ===" >> "${log_file}"
    rm -f "${sv_prep_ref_bam}"
    
    # Gridss 실행
    if [[ ! -f "${gridss_raw_vcf}" ]]; then
        log_info "[${sample_id}] Running Gridss"
        echo "=== GRIDSS START: $(date) ===" >> "${log_file}"
        
        "${GRIDSS_SCRIPT}" --jar "${GRIDSS_JAR}" \
            --steps all \
            --workingdir "${OUTPUT_DIR}" \
            --reference "${REF_GENOME}" \
            --blacklist "${GRIDSS_BLACKLIST_BED}" \
            --configuration "${GRIDSS_CONFIG}" \
            --labels "${sample_id}" \
            --bams "${reference_bam}" \
            --filtered_bams "${sv_prep_ref_sorted_bam}" \
            --output "${gridss_raw_vcf}" \
            --jvmheap "${sample_memory}G" \
            --threads "${GRIDSS_THREADS}" >> "${log_file}" 2>&1
        
        local gridss_exit_code=$?
        echo "=== GRIDSS END: $(date), Exit Code: ${gridss_exit_code} ===" >> "${log_file}"
        
        if [[ ! -f "${gridss_raw_vcf}" ]]; then
            error_exit "[${sample_id}] Gridss failed. Check log: ${log_file}"
        fi
    else
        log_info "[${sample_id}] Skipping Gridss, output exists"
    fi
    
    # SvPrep Reference Depth Annotator 실행
    log_info "[${sample_id}] Running SvPrep reference depth annotation"
    echo "=== DepthAnnotator START: $(date) ===" >> "${log_file}"
    
    java -Xmx${sample_memory}G -cp "${SV_PREP_JAR}" com.hartwig.hmftools.svprep.depth.DepthAnnotator \
        -input_vcf "${gridss_raw_vcf}" \
        -output_vcf "${gridss_vcf}" \
        -samples "${sample_id}" \
        -bam_files "${reference_bam}" \
        -ref_genome "${REF_GENOME}" \
        -ref_genome_version "${REF_GENOME_VERSION}" \
        -threads "${sample_threads}" >> "${log_file}" 2>&1
    
    local depth_exit_code=$?
    echo "=== DepthAnnotator END: $(date), Exit Code: ${depth_exit_code} ===" >> "${log_file}"
    
    if [[ ${depth_exit_code} -ne 0 ]]; then
        error_exit "[${sample_id}] DepthAnnotator failed. Check log: ${log_file}"
    fi
    
    log_info "[${sample_id}] GRIDSS completed successfully"
}

# GRIPSS 실행 함수
run_gripss() {
    local sample_id=$1
    local sample_memory=$2
    local log_file="${LOG_DIR}/gripss_${sample_id}_${TIMESTAMP}.log"
    
    log_info "[${sample_id}] Starting GRIPSS"
    
    local gridss_raw_vcf="${OUTPUT_DIR}/${sample_id}.gridss.raw.vcf.gz"
    local output_vcf="${OUTPUT_DIR}/${sample_id}.gripss.germline.vcf.gz"
    local filtered_vcf="${OUTPUT_DIR}/${sample_id}.gripss.filtered.germline.vcf.gz"
    
    if [[ -f "${filtered_vcf}" ]]; then
        log_info "[${sample_id}] Skipping GRIPSS, output exists"
    elif [[ ! -f "${gridss_raw_vcf}" ]]; then
        error_exit "[${sample_id}] Missing Gridss VCF: ${gridss_raw_vcf}"
    else
        log_info "[${sample_id}] Running GRIPSS"
        echo "=== GRIPSS START: $(date) ===" >> "${log_file}"
        
        java -Xmx${sample_memory}G -jar "${GRIPSS_JAR}" \
            -sample "${sample_id}" \
            -germline \
            -ref_genome "${REF_GENOME}" \
            -ref_genome_version "${REF_GENOME_VERSION}" \
            -known_hotspot_file "${SV_HOTSPOT}" \
            -pon_sv_file "${SV_PON_FILE}" \
            -pon_sgl_file "${SGL_PON_FILE}" \
            -repeat_mask_file "${REPEAT_MASK_FILE}" \
            -vcf "${gridss_raw_vcf}" \
            -output_dir "${OUTPUT_DIR}" \
            -output_id germline >> "${log_file}" 2>&1
        
        local exit_code=$?
        echo "=== GRIPSS END: $(date), Exit Code: ${exit_code} ===" >> "${log_file}"
        
        if [[ ${exit_code} -eq 0 ]]; then
            log_info "[${sample_id}] GRIPSS completed successfully"
        else
            error_exit "[${sample_id}] GRIPSS failed. Check log: ${log_file}"
        fi
    fi
}

# PURPLE 실행 함수
run_purple() {
    local sample_id=$1
    local sample_memory=$2
    local sample_threads=$3
    local log_file="${LOG_DIR}/purple_${sample_id}_${TIMESTAMP}.log"
    
    log_info "[${sample_id}] Starting PURPLE"
    
    local pave_vcf="${OUTPUT_DIR}/${sample_id}.pave.germline.vcf.gz"
    local gripss_vcf="${OUTPUT_DIR}/${sample_id}.gripss.filtered.germline.vcf.gz"
    
    # 입력 파일 확인
    if [[ ! -f "${pave_vcf}" ]]; then
        error_exit "[${sample_id}] Missing PAVE VCF: ${pave_vcf}"
    fi
    if [[ ! -f "${gripss_vcf}" ]]; then
        error_exit "[${sample_id}] Missing GRIPSS VCF: ${gripss_vcf}"
    fi
    
    log_info "[${sample_id}] Running PURPLE"
    echo "=== PURPLE START: $(date) ===" >> "${log_file}"
    
    java -Xmx${sample_memory}G \
        -jar "${PURPLE_JAR}" \
        -reference "${sample_id}" \
        -amber "${OUTPUT_DIR}" \
        -cobalt "${OUTPUT_DIR}" \
        -gc_profile "${GC_PROFILE}" \
        -ref_genome_version "${REF_GENOME_VERSION}" \
        -ref_genome "${REF_GENOME}" \
        -ensembl_data_dir "${ENSEMBL_DIR}" \
        -germline_hotspots "${GERMLINE_HOTSPOTS}" \
        -germline_del_freq_file "${GERMLINE_DEL_FREQ_FILE}" \
        -germline_vcf "${pave_vcf}" \
        -germline_sv_vcf "${gripss_vcf}" \
        -driver_gene_panel "${DRIVER_GENE_PANEL}" \
        -output_dir "${OUTPUT_DIR}" \
        -threads "${sample_threads}" >> "${log_file}" 2>&1
    
    local exit_code=$?
    echo "=== PURPLE END: $(date), Exit Code: ${exit_code} ===" >> "${log_file}"
    
    if [[ ${exit_code} -eq 0 ]]; then
        log_info "[${sample_id}] PURPLE completed successfully"
    else
        error_exit "[${sample_id}] PURPLE failed. Check log: ${log_file}"
    fi
}

# LINX 실행 함수
run_linx() {
    local sample_id=$1
    local sample_memory=$2
    local log_file="${LOG_DIR}/linx_${sample_id}_${TIMESTAMP}.log"
    
    log_info "[${sample_id}] Starting LINX"
    
    local purple_sv_vcf="${OUTPUT_DIR}/${sample_id}.purple.sv.germline.vcf.gz"
    
    if [[ ! -f "${purple_sv_vcf}" ]]; then
        error_exit "[${sample_id}] Missing PURPLE SV VCF: ${purple_sv_vcf}"
    fi
    
    log_info "[${sample_id}] Running LINX"
    echo "=== LINX START: $(date) ===" >> "${log_file}"
    
    java -Xmx${sample_memory}G \
        -jar "${LINX_JAR}" \
        -sample "${sample_id}" \
        -germline \
        -ref_genome_version "${REF_GENOME_VERSION}" \
        -sv_vcf "${purple_sv_vcf}" \
        -output_dir "${OUTPUT_DIR}" \
        -ensembl_data_dir "${ENSEMBL_DIR}" \
        -driver_gene_panel "${DRIVER_GENE_PANEL}" >> "${log_file}" 2>&1
    
    local exit_code=$?
    echo "=== LINX END: $(date), Exit Code: ${exit_code} ===" >> "${log_file}"
    
    if [[ ${exit_code} -eq 0 ]]; then
        log_info "[${sample_id}] LINX completed successfully"
    else
        error_exit "[${sample_id}] LINX failed. Check log: ${log_file}"
    fi
}

# 단일 샘플 전체 파이프라인 실행 (각 샘플의 main 함수)
run_sample_main() {
    local sample_id=$1
    local sample_memory=$2
    local sample_threads=$3
    
    log_info "=== Starting complete pipeline for sample: ${sample_id} ==="
    log_info "[${sample_id}] Memory: ${sample_memory}G, Threads: ${sample_threads}"
    
    # Stage 1: 병렬 실행 (SAGE_PAVE, AMBER_COBALT, GRIDSS)
    log_info "[${sample_id}] === Stage 1: SAGE_PAVE, AMBER_COBALT, GRIDSS (parallel) ==="
    
    run_sage_pave "${sample_id}" "${sample_memory}" "${sample_threads}" &
    local sage_pave_pid=$!
    
    run_amber_cobalt "${sample_id}" "${sample_memory}" "${sample_threads}" &
    local amber_cobalt_pid=$!
    
    run_gridss "${sample_id}" "${sample_memory}" "${sample_threads}" &
    local gridss_pid=$!
    
    # Stage 1 완료 대기
    wait ${sage_pave_pid} ${amber_cobalt_pid} ${gridss_pid}
    log_info "[${sample_id}] === Stage 1 completed ==="
    
    # Stage 2: 순차 실행 (GRIPSS, PURPLE)
    log_info "[${sample_id}] === Stage 2: GRIPSS, PURPLE (sequential) ==="
    
    run_gripss "${sample_id}" "${sample_memory}"
    run_purple "${sample_id}" "${sample_memory}" "${sample_threads}"
    
    log_info "[${sample_id}] === Stage 2 completed ==="
    
    # Stage 3: LINX 실행
    log_info "[${sample_id}] === Stage 3: LINX ==="
    
    run_linx "${sample_id}" "${sample_memory}"
    
    log_info "[${sample_id}] === Stage 3 completed ==="
    
    # 임시 파일 정리
    find "${OUTPUT_DIR}" -name "tmp_${sample_id}*" -type d -exec rm -rf {} + 2>/dev/null || true
    
    log_info "=== Sample ${sample_id} pipeline completed successfully ==="
}

# Dry run 실행
run_dry_run() {
    local sample_ids=("$@")
    
    echo "=== DRY RUN MODE ==="
    echo "Would process ${#sample_ids[@]} samples with complete pipeline each:"
    echo "Memory per sample: ${SAMPLE_MEMORY}G"
    echo "Threads per sample: ${SAMPLE_THREADS}"
    echo "Parallel samples: ${PARALLEL_SAMPLES}"
    echo ""
    
    for sample_id in "${sample_ids[@]}"; do
        echo "Sample: ${sample_id}"
        echo "  Complete pipeline: Stage 1 → Stage 2 → Stage 3"
        echo "  Stage 1 (parallel within sample):"
        echo "    - SAGE_PAVE: ${BAM_DIR}/${sample_id}_markdup.bam → ${OUTPUT_DIR}/${sample_id}.pave.germline.vcf.gz"
        echo "    - AMBER_COBALT: ${BAM_DIR}/${sample_id}_markdup.bam → ${OUTPUT_DIR}/${sample_id}.amber.baf.tsv.gz"
        echo "    - GRIDSS: ${BAM_DIR}/${sample_id}_markdup.bam → ${OUTPUT_DIR}/${sample_id}_gridss.vcf.gz"
        echo "  Stage 2 (sequential within sample):"
        echo "    - GRIPSS: ${OUTPUT_DIR}/${sample_id}.gridss.raw.vcf.gz → ${OUTPUT_DIR}/${sample_id}.gripss.filtered.germline.vcf.gz"
        echo "    - PURPLE: Multiple inputs → ${OUTPUT_DIR}/${sample_id}.purple.sv.germline.vcf.gz"
        echo "  Stage 3:"
        echo "    - LINX: ${OUTPUT_DIR}/${sample_id}.purple.sv.germline.vcf.gz → ${OUTPUT_DIR}/${sample_id}.linx.*.tsv"
        echo ""
    done
}

# 메인 함수
main() {
    # 기본 설정
    local requested_samples=()
    local run_all=true
    local dry_run=false
    
    # 기본값 설정
    PARALLEL_SAMPLES=4
    SAMPLE_MEMORY=${MAX_MEMORY}
    SAMPLE_THREADS=${THREADS}
    
    # 명령행 인수 파싱
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_usage
                exit 0
                ;;
            -l|--list)
                list_samples
                exit 0
                ;;
            -a|--all)
                run_all=true
                shift
                ;;
            --parallel)
                PARALLEL_SAMPLES="$2"
                shift 2
                ;;
            --max-memory)
                SAMPLE_MEMORY="$2"
                shift 2
                ;;
            --threads)
                SAMPLE_THREADS="$2"
                shift 2
                ;;
            --dry-run)
                dry_run=true
                shift
                ;;
            -*)
                log_error "Unknown option: $1"
                show_usage
                exit 1
                ;;
            *)
                requested_samples+=("$1")
                run_all=false
                shift
                ;;
        esac
    done
    
    log_info "=== Hereditary Cancer Germline Pipeline v4 Started ==="
    log_info "Each sample runs complete pipeline independently"
    log_info "Timestamp: ${TIMESTAMP}"
    log_info "Main Log: ${MAIN_LOG}"
    log_info "Parallel samples: ${PARALLEL_SAMPLES}"
    log_info "Memory per sample: ${SAMPLE_MEMORY}G"
    log_info "Threads per sample: ${SAMPLE_THREADS}"
    
    # 사전 확인
    check_bam_files
    
    # 처리할 샘플 결정
    local final_samples=()
    if [[ "${run_all}" == "true" ]]; then
        get_available_samples
        final_samples=("${available_samples[@]}")
    else
        validate_samples "${requested_samples[@]}"
        final_samples=("${requested_samples[@]}")
    fi
    
    log_info "Processing ${#final_samples[@]} samples: ${final_samples[*]}"
    
    # Dry run 모드
    if [[ "${dry_run}" == "true" ]]; then
        run_dry_run "${final_samples[@]}"
        exit 0
    fi
    
    # 각 샘플별로 독립적인 파이프라인 실행
    log_info "Starting ${#final_samples[@]} independent sample pipelines..."
    
    local active_jobs=0
    for sample_id in "${final_samples[@]}"; do
        # 병렬 작업 수 제한
        while [[ ${active_jobs} -ge ${PARALLEL_SAMPLES} ]]; do
            sleep 10
            active_jobs=$(jobs -r | wc -l)
        done
        
        # 각 샘플의 전체 파이프라인을 백그라운드에서 실행
        run_sample_main "${sample_id}" "${SAMPLE_MEMORY}" "${SAMPLE_THREADS}" &
        ((active_jobs++))
        
        log_info "Started pipeline for sample ${sample_id} (${active_jobs} active jobs)"
    done
    
    # 모든 샘플 파이프라인 완료 대기
    log_info "Waiting for all sample pipelines to complete..."
    wait
    
    log_info "=== All sample pipelines completed ==="
    
    # 전역 정리
    log_info "Cleaning up temporary files..."
    find "${OUTPUT_DIR}" -name "tmp_*" -type d -exec rm -rf {} + 2>/dev/null || true
    
    # 완료 요약
    log_info "=== Hereditary Cancer Germline Pipeline v4 Completed Successfully ==="
    log_info "Processed ${#final_samples[@]} samples: ${final_samples[*]}"
    log_info "Results are available in: ${OUTPUT_DIR}"
    log_info "Logs are available in: ${LOG_DIR}"
    log_info "Main log: ${MAIN_LOG}"
    
    # 결과 요약 생성
    generate_summary "${final_samples[@]}"
}

# 결과 요약 생성
generate_summary() {
    local processed_samples=("$@")
    local summary_log="${LOG_DIR}/pipeline_summary_${TIMESTAMP}.log"
    
    echo "=== Pipeline Execution Summary ===" > "${summary_log}"
    echo "Execution Time: $(date)" >> "${summary_log}"
    echo "Samples Processed: ${#processed_samples[@]}" >> "${summary_log}"
    echo "Sample IDs: ${processed_samples[*]}" >> "${summary_log}"
    echo "Parallel samples: ${PARALLEL_SAMPLES}" >> "${summary_log}"
    echo "Memory per sample: ${SAMPLE_MEMORY}G" >> "${summary_log}"
    echo "Threads per sample: ${SAMPLE_THREADS}" >> "${summary_log}"
    echo "" >> "${summary_log}"
    
    echo "=== Log Files ===" >> "${summary_log}"
    echo "Main Log: ${MAIN_LOG}" >> "${summary_log}"
    echo "Individual Tool Logs:" >> "${summary_log}"
    find "${LOG_DIR}" -name "*_${TIMESTAMP}.log" -type f | sort >> "${summary_log}"
    echo "" >> "${summary_log}"
    
    echo "=== Output Files ===" >> "${summary_log}"
    for sample_id in "${processed_samples[@]}"; do
        echo "Sample: ${sample_id}" >> "${summary_log}"
        find "${OUTPUT_DIR}" -name "${sample_id}*" -type f | sort >> "${summary_log}"
        echo "" >> "${summary_log}"
    done
    
    log_info "Pipeline summary generated: ${summary_log}"
}

# 스크립트 실행
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
