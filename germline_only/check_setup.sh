#!/bin/bash

# =============================================================================
# Pipeline Setup Verification Script
# =============================================================================

# 스크립트 디렉토리 경로
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 설정 파일 로드
source "${SCRIPT_DIR}/pipeline_config.sh"

# 색상 정의
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 로그 함수
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 카운터
error_count=0
warn_count=0

# 디렉토리 존재 확인
check_directory() {
    local dir_path=$1
    local dir_name=$2
    
    if [[ -d "${dir_path}" ]]; then
        log_info "${dir_name} 디렉토리 존재: ${dir_path}"
    else
        log_error "${dir_name} 디렉토리 없음: ${dir_path}"
        ((error_count++))
    fi
}

# 파일 존재 확인
check_file() {
    local file_path=$1
    local file_name=$2
    local required=$3
    
    if [[ -f "${file_path}" ]]; then
        log_info "${file_name} 파일 존재: ${file_path}"
    else
        if [[ "${required}" == "required" ]]; then
            log_error "${file_name} 파일 없음: ${file_path}"
            ((error_count++))
        else
            log_warn "${file_name} 파일 없음: ${file_path}"
            ((warn_count++))
        fi
    fi
}

# JAR 파일 실행 가능 확인
check_jar() {
    local jar_path=$1
    local jar_name=$2
    
    if [[ -f "${jar_path}" ]]; then
        if java -jar "${jar_path}" --help >/dev/null 2>&1 || java -jar "${jar_path}" -h >/dev/null 2>&1; then
            log_info "${jar_name} JAR 파일 실행 가능: ${jar_path}"
        else
            log_warn "${jar_name} JAR 파일 실행 확인 실패: ${jar_path}"
            ((warn_count++))
        fi
    else
        log_error "${jar_name} JAR 파일 없음: ${jar_path}"
        ((error_count++))
    fi
}

# BAM 파일 확인
check_bam_files() {
    log_info "BAM 파일 확인 중..."
    
    if [[ ! -d "${BAM_DIR}" ]]; then
        log_error "BAM 디렉토리 없음: ${BAM_DIR}"
        ((error_count++))
        return
    fi
    
    bam_files=($(find "${BAM_DIR}" -name "*.bam" 2>/dev/null))
    
    if [[ ${#bam_files[@]} -eq 0 ]]; then
        log_error "BAM 파일이 없습니다: ${BAM_DIR}"
        ((error_count++))
    else
        log_info "발견된 BAM 파일 수: ${#bam_files[@]}"
        
        # 각 BAM 파일의 인덱스 확인
        for bam_file in "${bam_files[@]}"; do
            sample_id=$(basename "${bam_file}" .bam)
            bai_file="${bam_file}.bai"
            
            if [[ -f "${bai_file}" ]]; then
                log_info "BAM 인덱스 존재: ${sample_id}"
            else
                log_error "BAM 인덱스 없음: ${sample_id} (${bai_file})"
                ((error_count++))
            fi
        done
    fi
}

# 메모리 확인
check_memory() {
    log_info "메모리 설정 확인 중..."
    
    # 시스템 메모리 확인
    total_mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    total_mem_gb=$((total_mem_kb / 1024 / 1024))
    
    log_info "시스템 총 메모리: ${total_mem_gb}GB"
    log_info "설정된 Java 힙 메모리: ${MAX_MEMORY}GB"
    
    if [[ ${MAX_MEMORY} -gt $((total_mem_gb - 16)) ]]; then
        log_warn "Java 힙 메모리가 시스템 메모리에 비해 클 수 있습니다"
        ((warn_count++))
    fi
}

# CPU 확인
check_cpu() {
    log_info "CPU 설정 확인 중..."
    
    cpu_cores=$(nproc)
    log_info "시스템 CPU 코어 수: ${cpu_cores}"
    log_info "설정된 스레드 수: ${THREADS}"
    
    if [[ ${THREADS} -gt ${cpu_cores} ]]; then
        log_warn "설정된 스레드 수가 CPU 코어 수보다 많습니다"
        ((warn_count++))
    fi
}

# 디스크 공간 확인
check_disk_space() {
    log_info "디스크 공간 확인 중..."
    
    if [[ -d "${OUTPUT_DIR}" ]]; then
        available_space=$(df "${OUTPUT_DIR}" | awk 'NR==2 {print $4}')
        available_gb=$((available_space / 1024 / 1024))
        
        log_info "출력 디렉토리 사용 가능 공간: ${available_gb}GB"
        
        if [[ ${available_gb} -lt 100 ]]; then
            log_warn "출력 디렉토리의 사용 가능 공간이 부족할 수 있습니다 (< 100GB)"
            ((warn_count++))
        fi
    fi
}

# 메인 검증 함수
main() {
    echo "=============================================="
    echo "Pipeline Setup Verification"
    echo "=============================================="
    
    log_info "설정 파일 로드됨: ${SCRIPT_DIR}/pipeline_config.sh"
    
    # 기본 디렉토리 확인
    log_info "기본 디렉토리 확인 중..."
    check_directory "${BASE_DIR}" "Base"
    check_directory "${BAM_DIR}" "BAM"
    check_directory "${OUTPUT_DIR}" "Output" 
    check_directory "${TOOLS_DIR}" "Tools"
    check_directory "${RESOURCE_DIR}" "Resource"
    check_directory "${ENSEMBL_DIR}" "Ensembl"
    
    # JAR 파일 확인
    log_info "JAR 파일 확인 중..."
    check_jar "${SAGE_JAR}" "SAGE"
    check_jar "${PAVE_JAR}" "PAVE"
    check_jar "${AMBER_JAR}" "AMBER"
    check_jar "${COBALT_JAR}" "COBALT"
    check_jar "${SV_PREP_JAR}" "SV-PREP"
    check_jar "${GRIDSS_JAR}" "GRIDSS"
    check_jar "${PURPLE_JAR}" "PURPLE"
    check_jar "${GRIPSS_JAR}" "GRIPSS"
    check_jar "${LINX_JAR}" "LINX"
    
    # 참조 파일 확인
    log_info "참조 파일 확인 중..."
    check_file "${REF_GENOME}" "Reference Genome" "required"
    check_file "${HOTSPOTS}" "Hotspots" "required"
    check_file "${PANEL_BED}" "Panel BED" "optional"
    check_file "${COVERAGE_BED}" "Coverage BED" "optional"
    check_file "${HIGH_CONFIDENCE_BED}" "High Confidence BED" "required"
    check_file "${MAPPABILITY_BED}" "Mappability BED" "required"
    check_file "${CLINVAR_VCF}" "ClinVar VCF" "required"
    check_file "${BLACKLIST_BED}" "Blacklist BED" "required"
    check_file "${BLACKLIST_VCF}" "Blacklist VCF" "required"
    check_file "${LOCI_FILE}" "Loci File" "required"
    check_file "${GC_PROFILE}" "GC Profile" "required"
    check_file "${DRIVER_GENE_PANEL}" "Driver Gene Panel" "required"
    check_file "${SV_BLACKLIST_BED}" "SV Blacklist BED" "required"
    check_file "${KNOWN_FUSION_BED}" "Known Fusion BED" "required"
    check_file "${GRIDSS_BLACKLIST_BED}" "GRIDSS Blacklist BED" "required"
    check_file "${GRIDSS_CONFIG}" "GRIDSS Config" "required"
    check_file "${GERMLINE_HOTSPOTS}" "Germline Hotspots" "required"
    check_file "${GERMLINE_DEL_FREQ_FILE}" "Germline Del Freq" "required"
    check_file "${SV_HOTSPOT}" "SV Hotspot" "required"
    check_file "${SV_PON_FILE}" "SV PON" "required"
    check_file "${SGL_PON_FILE}" "SGL PON" "required"
    check_file "${REPEAT_MASK_FILE}" "Repeat Mask" "required"
    
    # GRIDSS 스크립트 확인
    check_file "${GRIDSS_SCRIPT}" "GRIDSS Script" "required"
    
    # BAM 파일 확인
    check_bam_files
    
    # 시스템 리소스 확인
    check_memory
    check_cpu
    check_disk_space
    
    # 요약
    echo "=============================================="
    echo "검증 완료"
    echo "=============================================="
    
    if [[ ${error_count} -eq 0 ]] && [[ ${warn_count} -eq 0 ]]; then
        log_info "모든 검사 통과! 파이프라인 실행 준비가 완료되었습니다."
    elif [[ ${error_count} -eq 0 ]]; then
        log_warn "경고 ${warn_count}개가 있지만 파이프라인 실행이 가능합니다."
    else
        log_error "오류 ${error_count}개, 경고 ${warn_count}개가 발견되었습니다."
        log_error "오류를 수정한 후 파이프라인을 실행하세요."
        exit 1
    fi
    
    echo ""
    log_info "파이프라인 실행 명령:"
    echo "./run_germline_pipeline.sh"
}

# 스크립트 실행
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
