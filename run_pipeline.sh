#!/bin/bash

# =============================================================================
# Hereditary Cancer Germline Analysis Pipeline - Nextflow Runner
# =============================================================================

set -e

# Default parameters
PROFILE="standard"
RESUME=""
WORK_DIR="work"
CONFIG_FILE="nextflow.config"

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -p, --profile PROFILE     Execution profile (standard, cluster, docker, singularity) [default: standard]
  -r, --resume             Resume previous run
  -w, --work-dir DIR       Work directory [default: work]
  -c, --config FILE        Configuration file [default: nextflow.config]
  --bam-dir DIR            BAM files directory
  --output-dir DIR         Output directory
  --max-cpus NUM           Maximum CPUs to use
  --max-memory SIZE        Maximum memory to use (e.g., 512.GB)
  --threads NUM            Threads per process
  --memory SIZE            Memory per process (e.g., 64.GB)
  --max-forks NUM          Maximum parallel samples (default: 8)
  --redux-jar PATH         Path to REDUX JAR file
  --samtools PATH          Path to samtools executable
  -h, --help               Show this help message

EXAMPLES:
  # Basic run
  $0

  # Run with specific resources
  $0 --max-cpus 32 --max-memory 256.GB --threads 8 --memory 32.GB

  # Resume previous run
  $0 --resume

  # Run with custom directories
  $0 --bam-dir /path/to/bams --output-dir /path/to/output

  # Run with Docker
  $0 --profile docker

  # Run on cluster
  $0 --profile cluster

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -r|--resume)
            RESUME="-resume"
            shift
            ;;
        -w|--work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --bam-dir)
            BAM_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --max-cpus)
            MAX_CPUS="$2"
            shift 2
            ;;
        --max-memory)
            MAX_MEMORY="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --memory)
            MEMORY="$2"
            shift 2
            ;;
        --max-forks)
            MAX_FORKS="$2"
            shift 2
            ;;
        --redux-jar)
            REDUX_JAR="$2"
            shift 2
            ;;
        --samtools)
            SAMTOOLS="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Build Nextflow command
NF_CMD="nextflow run main.nf"
NF_CMD="$NF_CMD -profile $PROFILE"
NF_CMD="$NF_CMD -work-dir $WORK_DIR"
NF_CMD="$NF_CMD -c $CONFIG_FILE"

if [[ -n "$RESUME" ]]; then
    NF_CMD="$NF_CMD $RESUME"
fi

if [[ -n "$BAM_DIR" ]]; then
    NF_CMD="$NF_CMD --bam_dir $BAM_DIR"
fi

if [[ -n "$OUTPUT_DIR" ]]; then
    NF_CMD="$NF_CMD --output_dir $OUTPUT_DIR"
fi

if [[ -n "$MAX_CPUS" ]]; then
    NF_CMD="$NF_CMD --max_cpus $MAX_CPUS"
fi

if [[ -n "$MAX_MEMORY" ]]; then
    NF_CMD="$NF_CMD --max_memory $MAX_MEMORY"
fi

if [[ -n "$THREADS" ]]; then
    NF_CMD="$NF_CMD --threads $THREADS"
fi

if [[ -n "$MEMORY" ]]; then
    NF_CMD="$NF_CMD --memory $MEMORY"
fi

if [[ -n "$MAX_FORKS" ]]; then
    NF_CMD="$NF_CMD --max_forks $MAX_FORKS"
fi

if [[ -n "$REDUX_JAR" ]]; then
    NF_CMD="$NF_CMD --redux_jar $REDUX_JAR"
fi

if [[ -n "$SAMTOOLS" ]]; then
    NF_CMD="$NF_CMD --samtools $SAMTOOLS"
fi

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: Nextflow is not installed or not in PATH"
    echo "Please install Nextflow: curl -s https://get.nextflow.io | bash"
    exit 1
fi

# Check if main.nf exists
if [[ ! -f "main.nf" ]]; then
    echo "ERROR: main.nf not found in current directory"
    exit 1
fi

# Create output directory if it doesn't exist
if [[ -n "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
elif grep -q "output_dir" nextflow.config; then
    OUTPUT_DIR=$(grep "output_dir" nextflow.config | head -1 | sed 's/.*= *"\([^"]*\)".*/\1/')
    mkdir -p "$OUTPUT_DIR"
fi

echo "=== Hereditary Cancer Germline Pipeline - Nextflow ==="
echo "Profile: $PROFILE"
echo "Work directory: $WORK_DIR"
echo "Config file: $CONFIG_FILE"
echo "Command: $NF_CMD"
echo "=================================================="

# Run the pipeline
eval $NF_CMD

echo "=================================================="
echo "Pipeline completed!"
if [[ -n "$OUTPUT_DIR" ]]; then
    echo "Results are available in: $OUTPUT_DIR"
    echo "Pipeline reports are available in: $OUTPUT_DIR/pipeline_info/"
fi
echo "=================================================="
