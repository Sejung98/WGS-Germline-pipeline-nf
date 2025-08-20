#!/bin/bash

# 기본 설정
bam_dir="/home/ricky8419/09_Hereditary_cancer/02_bam/after_markdup"  # .bam 파일이 있는 디렉토리
output_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results/"  # 출력 디렉토리
amber_jar="/home/ricky8419/HMF_pipeline/tool/amber.jar"  # amber.jar 위치
cobalt_jar="/home/ricky8419/HMF_pipeline/tool/cobalt.jar"  # cobalt.jar 위치
loci_file="/home/ricky8419/HMF_pipeline/resource/38/copy_number/GermlineHetPon.38.vcf.gz"  # loci 파일 경로
gc_profile="/home/ricky8419/HMF_pipeline/resource/38/copy_number/GC_profile.1000bp.38.cnp"  # GC profile 파일 경로
threads=32  # 스레드 수
memory="256G"  # 메모리 설정

# .bam 파일 처리
for bam_file in "${bam_dir}"/*.bam; do
  # 샘플 이름 추출 (파일명에서 .bam 제거)
  sample_id=$(basename "${bam_file}" .bam)

  # 참조 샘플 ID 설정
  reference_id="${sample_id}"
  reference_bam="${bam_dir}/${reference_id}.bam"

  # Amber 명령 실행
  echo "Processing sample: ${sample_id}"
  java -Xmx${memory} \
    -jar "${amber_jar}" \
    -reference "${reference_id}" \
    -reference_bam "${reference_bam}" \
    -ref_genome_version 38 \
    -output_dir "${output_dir}" \
    -threads "${threads}" \
    -loci "${loci_file}"

  # 실행 결과 확인
  if [ $? -ne 0 ]; then
    echo "Error processing sample: ${sample_id}"
    exit 1
  fi

 # Cobalt 명령 실행
  echo "Running Cobalt for sample: ${sample_id}"

  java -Xmx${memory} \
    -jar "${cobalt_jar}" \
    -reference "${reference_id}" \
    -reference_bam "${reference_bam}" \
    -output_dir "${output_dir}" \
    -threads "${threads}" \
    -gc_profile "${gc_profile}"

  # Cobalt 실행 결과 확인
  if [ $? -ne 0 ]; then
    echo "Error processing Cobalt for sample: ${sample_id}"
    exit 1
  fi
done

echo "All samples processed successfully."

