# Hereditary Cancer Germline Pipeline

Integrated pipeline scripts for hereditary cancer germline analysis.

## File Structure

```
001_scripts/germline_only/
├── pipeline_config.sh              # Integrated configuration file
├── run_germline_pipeline_v4.sh     # Main execution script (independent pipeline per sample)
├── check_setup.sh                  # Configuration verification script
└── README.md                       # Usage guide
```

## Pipeline Stages

### Stage 1 (Parallel Execution)
- **SAGE + PAVE**: Variant calling and annotation
- **AMBER + COBALT**: Data preparation for copy number variant analysis
- **GRIDSS**: Structural variant calling

### Stage 2 (After Stage 1 Completion)
- **GRIPSS**: GRIDSS result filtering
- **PURPLE**: Integrated analysis of AMBER, COBALT, and PAVE results

### Stage 3 (After Stage 2 Completion)
- **LINX**: Structural variant analysis using PURPLE results

## Usage

### 1. Configuration Verification
First, verify that the paths in `pipeline_config.sh` are correct:

```bash
# Key configuration items
export BASE_DIR="/home/ricky8419/09_Hereditary_cancer"
export BAM_DIR="${BASE_DIR}/00_rawdata/bam_markdup"
export OUTPUT_DIR="${BASE_DIR}/Results_germlineMode"
export TOOLS_DIR="/home/ricky8419/HMF_pipeline/tool"
export RESOURCE_DIR="/home/ricky8419/HMF_pipeline/resource/38"
```

### 2. Pipeline Execution
```bash
cd /home/ricky8419/09_Hereditary_cancer/001_scripts/germline_only

# Process all samples (each sample runs independent complete pipeline)
./run_germline_pipeline_v4.sh

# Process specific samples only
./run_germline_pipeline_v4.sh sample1 sample2

# Adjust parallel sample count
./run_germline_pipeline_v4.sh --parallel 8

# Adjust resources per sample
./run_germline_pipeline_v4.sh --max-memory 128 --threads 16

# Check available sample list
./run_germline_pipeline_v4.sh --list

# Preview execution plan
./run_germline_pipeline_v4.sh --dry-run
```

**Key Features:**
- **Independent Pipeline per Sample**: Each sample runs the complete 1→2→3 stage pipeline independently
- **Parallel Sample Processing**: Process multiple samples simultaneously (default: 4)
- **Resource Optimization**: Memory/thread allocation per sample
- **Complete Independence**: One sample failure doesn't affect other samples
- **Flexible Configuration**: Adjustable parallel count, memory, and threads

### 3. Log Monitoring
Execution logs are stored in `${OUTPUT_DIR}/logs/` directory:

- `pipeline_main_YYYYMMDD_HHMMSS.log`: Main pipeline log
- `{tool}_{sample_id}_YYYYMMDD_HHMMSS.log`: Detailed logs for each tool (including start/end time and exit codes)
- `pipeline_summary_YYYYMMDD_HHMMSS.log`: Execution summary and output file list

## Input File Requirements

### BAM Files
- Location: `${BAM_DIR}/*_markdup.bam`
- Format: `{sample_id}_markdup.bam`
- Requirements: Index file (`.bam.bai`) required

### Reference Files
All reference files must be in the paths configured in `pipeline_config.sh`.

## Output Files

The following files are generated for each sample:

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

## System Requirements

- **CPU**: 64 cores (4 samples × 16 threads or 8 samples × 8 threads)
- **Memory**: 512GB (4 samples × 128GB or 8 samples × 64GB)
- **Disk**: Sufficient temporary storage space required (~100GB per sample recommended)
- **Java**: 8 or higher

**Resource Allocation Examples:**
```bash
# Process 4 samples simultaneously, 128GB per sample, 16 threads
./run_germline_pipeline_v4.sh --parallel 4 --max-memory 128 --threads 16

# Process 8 samples simultaneously, 64GB per sample, 8 threads
./run_germline_pipeline_v4.sh --parallel 8 --max-memory 64 --threads 8
```

## Configuration Customization

You can adjust the following items in `pipeline_config.sh`:

```bash
# System resources (for basic version)
export THREADS=64           # CPU thread count
export MAX_MEMORY=512       # Java heap memory (GB)
export GRIDSS_THREADS=8     # GRIDSS dedicated thread count

# Path settings
export BAM_DIR="..."        # BAM file directory
export OUTPUT_DIR="..."     # Output directory
export TOOLS_DIR="..."      # Tools directory
```

## SAGE Germline Mode Configuration

This pipeline is optimized to run SAGE in germline mode:

### Key Parameters
- `-germline`: Activate germline mode
- `-ref_sample_count 0`: Disable germline filter
- `-tumor` / `-tumor_bam`: Actually reference samples (labels are opposite in germline mode)
- **Whole Genome Analysis**: Remove panel_only to detect germline variants across the entire genome

### Used Resource Files
- **Hotspots**: `KnownHotspots.germline.38.vcf.gz` - ClinVar pathogenic/likely pathogenic variants
- **High Confidence**: GIAB high confidence regions
- **Blacklist**: Variants not to report (used in PAVE)

### PAVE Post-processing
PAVE is executed with germline settings to perform:
- ClinVar annotation addition
- Driver gene panel-based classification
- Blacklist variant filtering
- Mappability information addition

## Troubleshooting

### Common Issues

1. **Insufficient Memory**: Reduce `MAX_MEMORY` value or increase system memory.
2. **Insufficient Disk Space**: Ensure sufficient temporary storage space.
3. **File Permission Issues**: Verify read/write permissions for all input files and output directories.

### Log Review
Detailed error information for each stage can be found in log files in the `${LOG_DIR}` directory.

### Restart
The pipeline skips already completed tasks, so it's safe to restart after failure.

## Important Notes

- Do not modify files in the output directory during pipeline execution.
- Ensure sufficient disk space before execution.
- Adjust system resource settings according to hardware specifications.
