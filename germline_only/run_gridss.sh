#!/bin/bash

# 기본 설정
bam_dir="/home/ricky8419/09_Hereditary_cancer/02_bam/after_markdup"  # .bam 파일이 있는 디렉토리
output_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results"  # Gridss 출력 디렉토리
tools_dir="/home/ricky8419/HMF_pipeline/tool/"  # 도구 디렉토리
ref_genome_version="38"  # 참조 유전체 버전
ref_genome="/home/ricky8419/HMF_pipeline/resource/38/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"  # 참조 유전체 파일

sv_blacklist_bed="/home/ricky8419/HMF_pipeline/resource/38/sv/sv_prep_blacklist.38.bed"  # SV 블랙리스트 BED 파일
known_fusion_bed="/home/ricky8419/HMF_pipeline/resource/38/sv/known_fusions.38.bedpe"  # Known Fusion BED 파일
gridss_blacklist_bed="/home/ricky8419/HMF_pipeline/resource/38/sv/gridss_blacklist.38.bed"  # Gridss 블랙리스트 BED 파일
gridss_config="/home/ricky8419/HMF_pipeline/resource/38/sv/gridss.properties"  # Gridss 구성 파일

threads=32  # 스레드 수
max_memory=256  # Java 메모리 설정 (GB)

# 도구 파일 경로
sv_prep_jar="${tools_dir}/sv-prep.jar"
gridss_script="${tools_dir}/gridss.run.sh"
gridss_jar="${tools_dir}/gridss.jar"


# .bam 파일 처리
for bam_file in "${bam_dir}"/*.bam; do

  # 샘플 ID 추출
  reference_id=$(basename "${bam_file}" .bam)
  reference_bam="${bam_file}"
  sv_prep_ref_bam="${output_dir}/${reference_id}.sv_prep.bam"
  sv_prep_ref_sorted_bam="${output_dir}/${reference_id}.sv_prep.sorted.bam"
  gridss_raw_vcf="${output_dir}/${reference_id}.gridss.raw.vcf.gz"

  # SvPrep 실행
  sv_prep_write_types="JUNCTIONS;BAM;FRAGMENT_LENGTH_DIST"

  args="-sample ${reference_id} \
    -bam_file ${reference_bam} \
    -ref_genome ${ref_genome} \
    -ref_genome_version ${ref_genome_version} \
    -blacklist_bed ${sv_blacklist_bed} \
    -known_fusion_bed ${known_fusion_bed} \
    -write_types ${sv_prep_write_types} \
    -output_dir ${output_dir} \
    -threads ${threads}"

  echo "Running SvPrep on reference with args: ${args}"

  java -Xmx${max_memory}G -jar ${sv_prep_jar} ${args}

  # BAM 파일 정렬 및 인덱싱
  samtools sort -@ ${threads} -m 2G -T tmp -O bam ${sv_prep_ref_bam} -o ${sv_prep_ref_sorted_bam}
  samtools index -@ ${threads} ${sv_prep_ref_sorted_bam}

  rm ${sv_prep_ref_bam}

gridss_vcf="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results/${reference_id}_gridss.vcf.gz"  # Gridss VCF 파일 경로

  # Gridss 실행
  if [[ ! -f "${gridss_raw_vcf}" ]]; then
    args="--steps all \
      --workingdir ${output_dir} \
      --reference ${ref_genome} \
      --blacklist ${gridss_blacklist_bed} \
      --configuration ${gridss_config} \
      --labels ${reference_id} \
      --bams ${reference_bam} \
      --filtered_bams ${sv_prep_ref_sorted_bam} \
      --output ${gridss_raw_vcf} \
      --jvmheap ${max_memory}G \
      --threads 8"

    echo "Running Gridss with args: ${args}"
    ${gridss_script} --jar ${gridss_jar} ${args}

    if [[ ! -f "${gridss_raw_vcf}" ]]; then
      echo "Gridss failed - see error log:"
      cat ${output_dir}/gridss/gridss.full*.log
      exit 1
    fi
  else
    echo "Skipping Gridss process for ${reference_id}, since VCF ${gridss_raw_vcf} exists"
  fi

  # SvPrep로 Reference Depth Annotator 실행
  args="-input_vcf ${gridss_raw_vcf} \
    -output_vcf ${gridss_vcf} \
    -samples ${reference_id} \
    -bam_files ${reference_bam} \
    -ref_genome ${ref_genome} \
    -ref_genome_version ${ref_genome_version} \
    -threads ${threads}"

  echo "Running SvPrep reference depth annotation with args: ${args}"

  java -Xmx${max_memory}G -cp ${sv_prep_jar} com.hartwig.hmftools.svprep.depth.DepthAnnotator ${args}


# 모든 작업 후 삭제
if [[ -f "${sv_prep_ref_bam}" ]]; then
  rm ${sv_prep_ref_bam}
fi


done

echo "All samples processed successfully."
