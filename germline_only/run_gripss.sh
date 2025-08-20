#!/bin/bash

# 기본 설정
vcf_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results"  # Gridss VCF 파일이 있는 디렉토리
output_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results"  # Gripss 출력 디렉토리
tools_dir="/home/ricky8419/HMF_pipeline/tool/"  # 도구 디렉토리
ref_genome_version="38"  # 참조 유전체 버전
ref_genome="/home/ricky8419/HMF_pipeline/resource/38/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"  # 참조 유전체 파일

sv_hotspot="/home/ricky8419/HMF_pipeline/resource/38/sv/known_fusions.38.bedpe"  # SV 블랙리스트 BED 파일
sv_pon_file="/home/ricky8419/HMF_pipeline/resource/38/sv/sv_pon.38.bedpe.gz"
sgl_pon_file="/home/ricky8419/HMF_pipeline/resource/38/sv/sgl_pon.38.bed.gz"
repeat_mask_file="/home/ricky8419/HMF_pipeline/resource/38/sv/repeat_mask_data.38.fa.gz"

threads=8  # 스레드 수
max_memory=256  # Java 메모리 설정 (GB)

# 도구 파일 경로
gripss_jar="${tools_dir}/gripss.jar"


# .bam 파일 처리
for gridss_vcf in "${vcf_dir}"/*.gridss.raw.vcf.gz; do

  # 샘플 ID 추출
  reference_id=$(basename "${gridss_vcf}" .gridss.raw.vcf.gz)
  output_vcf="${output_dir}/${reference_id}.gripss.germline.vcf.gz"


  # Gripss 실행 스킵 조건

  if [[ -f "${output_vcf}" ]]; then
    echo "Skipping Gripss germline for ${reference_id}, since VCF ${output_vcf} exists"
    continue
  fi

  if [[ ! -f "${gridss_vcf}" ]]; then
    echo "Missing Gridss VCF, not running Gripss germline for ${reference_id}"
    continue
  fi

 # Gripss 실행 매개변수 구성
  args="-sample ${reference_id} \
    -germline \
    -ref_genome ${ref_genome} \
    -ref_genome_version ${ref_genome_version} \
    -known_hotspot_file ${sv_hotspot} \
    -pon_sv_file ${sv_pon_file} \
    -pon_sgl_file ${sgl_pon_file} \
    -repeat_mask_file ${repeat_mask_file} \
    -vcf ${gridss_vcf} \
    -output_dir ${output_dir} \
    -output_id germline"

  # Gripss 실행
  echo "Running Gripss germline for ${tumor_id} with args: ${args}"
  java -Xmx${max_memory}G -jar ${gripss_jar} ${args}

  # 실행 결과 확인
  if [[ $? -eq 0 ]]; then
    echo "Gripss germline completed successfully for ${tumor_id}. Output: ${output_vcf}"
  else
    echo "Error running Gripss germline for ${tumor_id}"
    exit 1
  fi
done

echo "All samples processed successfully."
