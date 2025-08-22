#!/bin/bash

# =============================================================================
# SAGE Tool Execution Script
# Germline variant calling with REDUX BAM
# =============================================================================

set -e

# 설정 파일 로드
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/pipeline_config.sh"

# 사용법 출력
show_usage() {
    cat << EOF
Usage: $0 [OPTIONS] <SAMPLE_ID>

SAGE 도구를 실행하여 germline variant calling을 수행합니다.
REDUX BAM 파일과 jitter 파라미터를 사용합니다.

OPTIONS:
  -h, --help              Show this help message
  -o, --output-dir DIR    Output directory (default: from config)
  -m, --memory MEM        Memory in GB (default: 64)
  -t, --threads NUM       Threads (default: 8)
  -b, --bam-file FILE     Input BAM file (default: REDUX BAM)
  --panel-only            Use panel-only mode (default: false)
  --dry-run               Show command without execution

EXAMPLES:
  $0 sample1                    # 기본 설정으로 실행 (REDUX BAM 사용)
  $0 sample1 -b custom.bam      # 커스텀 BAM 파일 사용
  $0 sample1 -m 128 -t 16      # 128GB 메모리, 16 스레드
  $0 sample1 --panel-only       # Panel-only 모드
  $0 sample1 --dry-run          # 실행할 명령어만 출력

OUTPUT FILES:
  - {SAMPLE_ID}.sage.germline.vcf.gz  # SAGE VCF 파일
EOF
}

# 기본값 설정
SAMPLE_ID=""
OUTPUT_DIR="${OUTPUT_DIR}"
MEMORY=64
THREADS=8
INPUT_BAM=""
PANEL_ONLY=false
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
        -b|--bam-file)
            INPUT_BAM="$2"
            shift 2
            ;;
        --panel-only)
            PANEL_ONLY=true
            shift
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

# BAM 파일 결정 (REDUX BAM 우선, 없으면 원본 BAM)
if [[ -z "${INPUT_BAM}" ]]; then
    REDUX_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.redux.bam"
    ORIGINAL_BAM="${BAM_DIR}/${SAMPLE_ID}_markdup.bam"
    
    if [[ -f "${REDUX_BAM}" ]]; then
        INPUT_BAM="${REDUX_BAM}"
        echo "Using REDUX BAM: ${INPUT_BAM}"
    elif [[ -f "${ORIGINAL_BAM}" ]]; then
        INPUT_BAM="${ORIGINAL_BAM}"
        echo "Using original BAM: ${INPUT_BAM}"
    else
        echo "ERROR: No BAM file found. Please provide with -b option or ensure REDUX BAM exists."
        exit 1
    fi
fi

# 입력 파일 확인
if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM file not found: ${INPUT_BAM}"
    exit 1
fi

# 출력 디렉토리 생성
mkdir -p "${OUTPUT_DIR}"

# 로그 파일 설정
LOG_FILE="${OUTPUT_DIR}/sage_${SAMPLE_ID}_$(date +%Y%m%d_%H%M%S).log"

# 출력 파일 경로
SAGE_VCF="${OUTPUT_DIR}/${SAMPLE_ID}.sage.germline.vcf.gz"

echo "=== SAGE Tool Execution ==="
echo "Sample ID: ${SAMPLE_ID}"
echo "Input BAM: ${INPUT_BAM}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Memory: ${MEMORY}GB"
echo "Threads: ${THREADS}"
echo "Panel only: ${PANEL_ONLY}"
echo "Log file: ${LOG_FILE}"
echo "================================"

# SAGE 명령어 구성
SAGE_CMD="java -Xmx${MEMORY}G -jar ${SAGE_JAR} \\
    -tumor \"${SAMPLE_ID}\" \\
    -tumor_bam \"${INPUT_BAM}\" \\
    -hotspots \"${HOTSPOTS}\" \\
    -high_confidence_bed \"${HIGH_CONFIDENCE_BED}\" \\
    -ref_genome \"${REF_GENOME}\" \\
    -ref_genome_version \"${REF_GENOME_VERSION}\" \\
    -ensembl_data_dir \"${ENSEMBL_DIR}\" \\
    -ref_sample_count 0 \\
    -germline"

# Panel-only 모드 설정
if [[ "${PANEL_ONLY}" == "true" ]]; then
    SAGE_CMD="${SAGE_CMD} \\
    -panel_bed \"${PANEL_BED}\" \\
    -panel_only"
fi

# Jitter 파라미터 디렉토리 추가 (REDUX BAM 사용 시)
if [[ "${INPUT_BAM}" == *".redux.bam" ]]; then
    SAGE_CMD="${SAGE_CMD} \\
    -jitter_param_dir \"${OUTPUT_DIR}\""
fi

# 나머지 파라미터 추가
SAGE_CMD="${SAGE_CMD} \\
    -threads ${THREADS} \\
    -output_vcf \"${SAGE_VCF}\""

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "=== DRY RUN MODE ==="
    echo "Command to execute:"
    echo "${SAGE_CMD}"
    exit 0
fi

# SAGE 실행
echo "Starting SAGE execution..."
echo "Command: ${SAGE_CMD}"
echo "Log file: ${LOG_FILE}"

# 실행 및 로깅
{
    echo "=== SAGE START: $(date) ==="
    echo "Sample: ${SAMPLE_ID}"
    echo "Input BAM: ${INPUT_BAM}"
    echo "Command: ${SAGE_CMD}"
    echo ""
    
    eval "${SAGE_CMD}"
    
    local exit_code=$?
    echo ""
    echo "=== SAGE END: $(date), Exit Code: ${exit_code} ==="
    
    if [[ ${exit_code} -eq 0 ]]; then
        echo "SAGE completed successfully"
        
        # 출력 파일 검증
        if [[ -f "${SAGE_VCF}" ]]; then
            echo "✓ SAGE VCF created: ${SAGE_VCF}"
            
            # VCF 인덱스 생성
            if [[ ! -f "${SAGE_VCF}.tbi" ]]; then
                echo "Creating VCF index..."
                "${TOOLS_DIR}/samtools" index "${SAGE_VCF}"
            fi
        else
            echo "✗ ERROR: SAGE VCF not created"
            exit 1
        fi
        
        echo "SAGE execution completed successfully!"
    else
        echo "SAGE failed with exit code ${exit_code}"
        exit ${exit_code}
    fi
} 2>&1 | tee "${LOG_FILE}"

echo ""
echo "SAGE execution completed. Check log file: ${LOG_FILE}"
echo "Output files:"
echo "  - ${SAGE_VCF}"
