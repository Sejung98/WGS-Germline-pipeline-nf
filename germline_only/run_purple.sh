#!/bin/bash

# 기본 설정
bam_dir="/home/ricky8419/09_Hereditary_cancer/02_bam/after_markdup"
output_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results"  # Purple 출력 디렉토리
tools_dir="/home/ricky8419/HMF_pipeline/tool/"  # 도구 디렉토리
ref_genome_version="38"  # 참조 유전체 버전
ref_genome="/home/ricky8419/HMF_pipeline/resource/38/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"  # 참조 유전체 파일

gc_profile="/home/ricky8419/HMF_pipeline/resource/38/copy_number/GC_profile.1000bp.38.cnp"  # GC Profile 파일
ensembl_dir="/home/ricky8419/HMF_pipeline/resource/38/common/ensembl_data"  # Ensembl 데이터 디렉토리
germline_hotspots="/home/ricky8419/HMF_pipeline/resource/38/variants/KnownHotspots.germline.38.vcf.gz"  # Germline Hotspots VCF
germline_del_freq_file="/home/ricky8419/HMF_pipeline/resource/38/copy_number/cohort_germline_del_freq.38.csv"  # Germline Deletion Frequency 파일
driver_gene_panel="/home/ricky8419/HMF_pipeline/resource/38/common/DriverGenePanel.38.tsv"  # Driver Gene Panel 파일

amber_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results"  # Amber 출력 디렉토리
cobalt_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results"  # Cobalt 출력 디렉토리

threads=32  # 스레드 수
max_memory=256  # Java 메모리 설정 (GB)

# Purple jar 파일 경로
purple_jar="${tools_dir}/purple.jar"

# 출력 디렉토리 생성
if [[ ! -d "${output_dir}" ]]; then
  mkdir -p ${output_dir}
fi

# .bam 파일 처리
for bam_file in "${bam_dir}"/*.bam; do
  # 샘플 ID 추출
  reference_id=$(basename "${bam_file}" .bam)


  # Purple 실행
  echo "Running Purple Germline-only for ${reference_id}..."

  java -Xmx${max_memory}G \
    -jar ${purple_jar} \
    -reference ${reference_id} \
    -amber ${amber_dir} \
    -cobalt ${cobalt_dir} \
    -gc_profile ${gc_profile} \
    -ref_genome_version ${ref_genome_version} \
    -ref_genome ${ref_genome} \
    -ensembl_data_dir ${ensembl_dir} \
    -germline_hotspots ${germline_hotspots} \
    -germline_del_freq_file ${germline_del_freq_file} \
    -germline_vcf ${output_dir}/${reference_id}.pave.germline.vcf.gz \
    -germline_sv_vcf ${output_dir}/${reference_id}.gripss.filtered.germline.vcf.gz \
    -driver_gene_panel ${driver_gene_panel} \
    -output_dir ${output_dir} \
    -threads 32

  # 실행 결과 확인
  if [[ $? -eq 0 ]]; then
    echo "Purple completed successfully for ${reference_id}. Output directory: ${output_dir}/${reference_id}"
  else
    echo "Error running Purple for ${reference_id}. Check logs for details."
    exit 1
  fi
done

echo "All samples processed successfully."

