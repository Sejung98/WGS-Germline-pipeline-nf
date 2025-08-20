#!/bin/bash

# 기본 설정
output_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results"  # LINX 출력 디렉토리
tools_dir="/home/ricky8419/HMF_pipeline/tool/"  # 도구 디렉토리
ref_genome_version="38"  # 참조 유전체 버전
ensembl_dir="/home/ricky8419/HMF_pipeline/resource/38/common/ensembl_data"  # Ensembl 데이터 경로
driver_gene_panel="/home/ricky8419/HMF_pipeline/resource/38/common/DriverGenePanel.38.tsv"  # Driver Gene Panel 경로

# LINX jar 파일 경로
linx_jar="${tools_dir}/linx.jar"

# SV VCF 파일 위치
purple_vcf_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results"  # Purple에서 생성된 VCF 디렉토리


# Purple VCF 처리
for sv_vcf in "${purple_vcf_dir}"/*.purple.sv.germline.vcf.gz; do
  # 샘플 ID 추출
  sample_id=$(basename "${sv_vcf}" .purple.sv.germline.vcf.gz)

  # LINX 실행
  echo "Running LINX for ${sample_id}..."

  java -Xmx256G \
    -jar ${linx_jar} \
    -sample ${sample_id} \
    -germline \
    -ref_genome_version ${ref_genome_version} \
    -sv_vcf ${sv_vcf} \
    -output_dir ${output_dir} \
    -ensembl_data_dir ${ensembl_dir} \
    -driver_gene_panel ${driver_gene_panel}

  # 실행 결과 확인
  if [[ $? -eq 0 ]]; then
    echo "LINX completed successfully for ${sample_id}. Output directory: ${output_dir}/${sample_id}"
  else
    echo "Error running LINX for ${sample_id}. Check logs for details."
    exit 1
  fi
done

echo "All samples processed successfully."
