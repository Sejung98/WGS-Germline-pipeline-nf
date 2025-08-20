# Hereditary Cancer Germline Pipeline

유전성 암 germline 분석을 위한 통합 파이프라인 스크립트입니다.

## 파일 구조

```
001_scripts/germline_only/
├── pipeline_config.sh              # 통합 설정 파일
├── run_germline_pipeline_v4.sh     # 메인 실행 스크립트 (샘플별 독립 파이프라인)
├── check_setup.sh                  # 설정 검증 스크립트
└── README.md                       # 사용 가이드
```

## 파이프라인 단계

### 1단계 (병렬 실행)
- **SAGE + PAVE**: 변이 호출 및 주석
- **AMBER + COBALT**: 복사수 변이 분석을 위한 데이터 준비
- **GRIDSS**: 구조적 변이 호출

### 2단계 (1단계 완료 후)
- **GRIPSS**: GRIDSS 결과 필터링
- **PURPLE**: AMBER, COBALT, PAVE 결과 통합 분석

### 3단계 (2단계 완료 후)
- **LINX**: PURPLE 결과로 구조적 변이 분석

## 사용법

### 1. 설정 확인
먼저 `pipeline_config.sh` 파일에서 경로들이 올바른지 확인하세요:

```bash
# 주요 설정 항목들
export BASE_DIR="/home/ricky8419/09_Hereditary_cancer"
export BAM_DIR="${BASE_DIR}/00_rawdata/bam_markdup"
export OUTPUT_DIR="${BASE_DIR}/Results_germlineMode"
export TOOLS_DIR="/home/ricky8419/HMF_pipeline/tool"
export RESOURCE_DIR="/home/ricky8419/HMF_pipeline/resource/38"
```

### 2. 파이프라인 실행
```bash
cd /home/ricky8419/09_Hereditary_cancer/001_scripts/germline_only

# 모든 샘플 처리 (각 샘플이 독립적인 전체 파이프라인 실행)
./run_germline_pipeline_v4.sh

# 특정 샘플만 처리
./run_germline_pipeline_v4.sh sample1 sample2

# 병렬 샘플 수 조정
./run_germline_pipeline_v4.sh --parallel 8

# 샘플당 리소스 조정
./run_germline_pipeline_v4.sh --max-memory 128 --threads 16

# 사용 가능한 샘플 목록 확인
./run_germline_pipeline_v4.sh --list

# 실행 계획 미리보기
./run_germline_pipeline_v4.sh --dry-run
```

**주요 특징:**
- **샘플별 독립 파이프라인**: 각 샘플이 1→2→3단계 전체 파이프라인을 독립적으로 실행
- **병렬 샘플 처리**: 여러 샘플을 동시에 처리 (기본 4개)
- **리소스 최적화**: 샘플당 메모리/스레드 할당
- **완전 독립성**: 한 샘플 실패가 다른 샘플에 영향 없음
- **유연한 설정**: 병렬 수, 메모리, 스레드 조정 가능

### 3. 로그 확인
실행 로그는 `${OUTPUT_DIR}/logs/` 디렉토리에 저장됩니다:

- `pipeline_main_YYYYMMDD_HHMMSS.log`: 메인 파이프라인 로그
- `{tool}_{sample_id}_YYYYMMDD_HHMMSS.log`: 각 도구별 상세 로그 (시작/종료 시간, 종료 코드 포함)
- `pipeline_summary_YYYYMMDD_HHMMSS.log`: 실행 요약 및 출력 파일 목록

## 입력 파일 요구사항

### BAM 파일
- 위치: `${BAM_DIR}/*_markdup.bam`
- 형식: `{sample_id}_markdup.bam`
- 요구사항: 인덱스 파일 (`.bam.bai`) 필요

### 참조 파일들
모든 참조 파일들은 `pipeline_config.sh`에서 설정된 경로에 있어야 합니다.

## 출력 파일

각 샘플에 대해 다음 파일들이 생성됩니다:

### SAGE & PAVE
- `{sample_id}.sage.germline.vcf.gz`
- `{sample_id}.pave.germline.vcf.gz`

### AMBER & COBALT
- `{sample_id}.amber.baf.tsv.gz`
- `{sample_id}.cobalt.ratio.tsv.gz`

### GRIDSS
- `{sample_id}.gridss.raw.vcf.gz`
- `{sample_id}_gridss.vcf.gz`

### GRIPSS
- `{sample_id}.gripss.germline.vcf.gz`
- `{sample_id}.gripss.filtered.germline.vcf.gz`

### PURPLE
- `{sample_id}.purple.sv.germline.vcf.gz`
- `{sample_id}.purple.cnv.germline.tsv`
- `{sample_id}.purple.driver.catalog.germline.tsv`

### LINX
- `{sample_id}.linx.germline.disruption.tsv`
- `{sample_id}.linx.germline.breakend.tsv`

## 시스템 요구사항

- **CPU**: 64 코어 (4개 샘플 × 16 스레드 또는 8개 샘플 × 8 스레드)
- **메모리**: 512GB (4개 샘플 × 128GB 또는 8개 샘플 × 64GB)
- **디스크**: 충분한 임시 저장 공간 필요 (샘플당 ~100GB 권장)
- **Java**: 8 이상

**리소스 할당 예시:**
```bash
# 4개 샘플 동시 처리, 샘플당 128GB, 16 스레드
./run_germline_pipeline_v4.sh --parallel 4 --max-memory 128 --threads 16

# 8개 샘플 동시 처리, 샘플당 64GB, 8 스레드  
./run_germline_pipeline_v4.sh --parallel 8 --max-memory 64 --threads 8
```

## 설정 커스터마이징

`pipeline_config.sh`에서 다음 항목들을 조정할 수 있습니다:

```bash
# 시스템 리소스 (기본 버전용)
export THREADS=64           # CPU 스레드 수
export MAX_MEMORY=512       # Java 힙 메모리 (GB)
export GRIDSS_THREADS=8     # GRIDSS 전용 스레드 수

# 경로 설정
export BAM_DIR="..."        # BAM 파일 디렉토리
export OUTPUT_DIR="..."     # 출력 디렉토리
export TOOLS_DIR="..."      # 도구 디렉토리
```

## SAGE Germline Mode 설정

이 파이프라인은 SAGE를 germline mode로 실행하도록 최적화되어 있습니다:

### 주요 파라미터
- `-germline`: germline mode 활성화
- `-ref_sample_count 0`: germline 필터 비활성화
- `-tumor` / `-tumor_bam`: 실제로는 reference 샘플 (germline mode에서는 라벨이 반대)
- **전체 게놈 분석**: panel_only 제거하여 전체 게놈에서 germline 변이 검출

### 사용된 리소스 파일
- **Hotspots**: `KnownHotspots.germline.38.vcf.gz` - ClinVar pathogenic/likely pathogenic 변이
- **High Confidence**: GIAB high confidence 영역
- **Blacklist**: 보고하지 않을 변이 목록 (PAVE에서 사용)

### PAVE 후처리
PAVE는 germline 설정으로 실행되어 다음을 수행합니다:
- ClinVar 주석 추가
- Driver gene panel 기반 분류
- Blacklist 변이 필터링
- Mappability 정보 추가

## 문제 해결

### 일반적인 문제들

1. **메모리 부족**: `MAX_MEMORY` 값을 줄이거나 시스템 메모리를 늘리세요.
2. **디스크 공간 부족**: 충분한 임시 저장 공간을 확보하세요.
3. **파일 권한 문제**: 모든 입력 파일과 출력 디렉토리에 대한 읽기/쓰기 권한을 확인하세요.

### 로그 확인
각 단계별 상세한 오류 정보는 `${LOG_DIR}` 디렉토리의 로그 파일에서 확인할 수 있습니다.

### 재시작
파이프라인은 이미 완료된 작업을 건너뛰므로, 실패 후 재시작해도 안전합니다.

## 주의사항

- 파이프라인 실행 중에는 출력 디렉토리의 파일을 수정하지 마세요.
- 충분한 디스크 공간을 확보한 후 실행하세요.
- 시스템 리소스 설정을 하드웨어 사양에 맞게 조정하세요.
