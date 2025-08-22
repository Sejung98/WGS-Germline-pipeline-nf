#!/bin/bash

# =============================================================================
# REDUX Tool Execution Script
# BAM 정제 및 microsatellite jitter 모델링
# =============================================================================

set -e

# 설정 파일 로드
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/pipeline_config.sh"

# 사용법 출력
show_usage() {
    cat << EOF
Usage: $0 [OPTIONS] <SAMPLE_ID>

REDUX 도구를 실행하여 BAM 파일을 정제하고 jitter 파라미터를 생성합니다.

OPTIONS:
  -h, --help              Show this help message
  -o, --output-dir DIR    Output directory (default: from config)
  -m, --memory MEM        Memory in GB (default: 64)
  -t, --threads NUM       Threads (default: 8)
  --dry-run               Show command without execution

EXAMPLES:
  $0 sample1                    # 기본 설정으로 실행
  $0 sample1 -o /custom/output  # 커스텀 출력 디렉토리
  $0 sample1 -m 128 -t 16      # 128GB 메모리, 16 스레드
  $0 sample1 --dry-run          # 실행할 명령어만 출력

OUTPUT FILES:
  - {SAMPLE_ID}.redux.bam           # 정제된 BAM 파일
  - {SAMPLE_ID}.jitter_params.tsv   # Jitter 파라미터
  - {SAMPLE_ID}.ms_table.tsv.gz     # Microsatellite 테이블
EOF
}

# 기본값 설정
SAMPLE_ID=""
OUTPUT_DIR="${OUTPUT_DIR}"
MEMORY=64
THREADS=8
DRY_RUN=false

# 명령행 인수 파싱
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -*)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
        *)
            if [[ -z "${SAMPLE_ID}" ]]; then
                SAMPLE_ID="$1"
            else
                echo "Multiple sample IDs provided. Only one allowed."
                exit 1
            fi
            shift
            ;;
    esac
done

# 샘플 ID 검증
if [[ -z "${SAMPLE_ID}" ]]; then
    echo "ERROR: Sample ID is required"
    show_usage
    exit 1
fi

# 입력 파일 확인
INPUT_BAM="${BAM_DIR}/${SAMPLE_ID}_markdup.bam"
if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM file not found: ${INPUT_BAM}"
    exit 1
fi

# 출력 디렉토리 생성
mkdir -p "${OUTPUT_DIR}"

# 로그 파일 설정
LOG_FILE="${OUTPUT_DIR}/redux_${SAMPLE_ID}_$(date +%Y%m%d_%H%M%S).log"

# 출력 파일 경로
REDUX_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.redux.bam"
JITTER_PARAMS="${OUTPUT_DIR}/${SAMPLE_ID}.jitter_params.tsv"
MS_TABLE="${OUTPUT_DIR}/${SAMPLE_ID}.ms_table.tsv.gz"

echo "=== REDUX Tool Execution ==="
echo "Sample ID: ${SAMPLE_ID}"
echo "Input BAM: ${INPUT_BAM}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Memory: ${MEMORY}GB"
echo "Threads: ${THREADS}"
echo "Log file: ${LOG_FILE}"
echo "================================"

# REDUX 명령어 구성
REDUX_CMD="java -Xmx${MEMORY}G -jar ${REDUX_JAR} \\
    -sample \"${SAMPLE_ID}\" \\
    -input_bam \"${INPUT_BAM}\" \\
    -output_bam \"${REDUX_BAM}\" \\
    -ref_genome \"${REF_GENOME}\" \\
    -ref_genome_version \"V${REF_GENOME_VERSION}\" \\
    -unmap_regions \"${RESOURCE_DIR}/common/unmap_regions.38.tsv\" \\
    -ref_genome_msi_file \"${RESOURCE_DIR}/common/msi_jitter_sites.38.tsv.gz\" \\
    -form_consensus \\
    -write_stats \\
    -bamtool \"${TOOLS_DIR}/samtools\" \\
    -output_dir \"${OUTPUT_DIR}\" \\
    -threads ${THREADS} \\
    -log_level INFO"

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "=== DRY RUN MODE ==="
    echo "Command to execute:"
    echo "${REDUX_CMD}"
    exit 0
fi

# REDUX 실행
echo "Starting REDUX execution..."
echo "Command: ${REDUX_CMD}"
echo "Log file: ${LOG_FILE}"

# 실행 및 로깅
{
    echo "=== REDUX START: $(date) ==="
    echo "Sample: ${SAMPLE_ID}"
    echo "Input BAM: ${INPUT_BAM}"
    echo "Command: ${REDUX_CMD}"
    echo ""
    
    eval "${REDUX_CMD}"
    
    local exit_code=$?
    echo ""
    echo "=== REDUX END: $(date), Exit Code: ${exit_code} ==="
    
    if [[ ${exit_code} -eq 0 ]]; then
        echo "REDUX completed successfully"
        
        # 출력 파일 검증
        if [[ -f "${REDUX_BAM}" ]]; then
            echo "✓ REDUX BAM created: ${REDUX_BAM}"
        else
            echo "✗ ERROR: REDUX BAM not created"
        fi
        
        if [[ -f "${JITTER_PARAMS}" ]]; then
            echo "✓ Jitter parameters created: ${JITTER_PARAMS}"
        else
            echo "✗ ERROR: Jitter parameters not created"
        fi
        
        if [[ -f "${MS_TABLE}" ]]; then
            echo "✓ MS table created: ${MS_TABLE}"
        else
            echo "✗ ERROR: MS table not created"
        fi
        
        # BAM 인덱스 생성
        echo "Creating BAM index..."
        "${TOOLS_DIR}/samtools" index "${REDUX_BAM}"
        
        echo "All REDUX outputs created successfully!"
    else
        echo "REDUX failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
} 2>&1 | tee "${LOG_FILE}"

echo ""
echo "REDUX execution completed. Check log file: ${LOG_FILE}"
echo "Output files:"
echo "  - ${REDUX_BAM}"
echo "  - ${JITTER_PARAMS}"
echo "  - ${MS_TABLE}"
