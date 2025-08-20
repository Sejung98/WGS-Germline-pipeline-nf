#!/bin/bash

# 기본 설정
bam_dir="/home/ricky8419/09_Hereditary_cancer/02_bam/after_markdup"  # .bam 파일이 있는 디렉토리
output_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results/"  # SAGE 출력 디렉토리
pave_output_dir="/home/ricky8419/09_Hereditary_cancer/03_Germ_only_results/"  # PAVE 출력 디렉토리

sage_jar="/home/ricky8419/HMF_pipeline/tool/sage.jar"  # SAGE jar 파일 경로
pave_jar="/home/ricky8419/HMF_pipeline/tool/pave.jar"  # PAVE jar 파일 경로

ref_genome_version="38"  # 참조 유전체 버전
ref_genome="/home/ricky8419/HMF_pipeline/resource/38/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"  # 참조 유전체 파일
ensembl_dir="/home/ricky8419/HMF_pipeline/resource/38/common/ensembl_data"  # Ensembl 데이터 경로
hotspots="/home/ricky8419/HMF_pipeline/resource/38/variants/KnownHotspots.germline.38.vcf.gz"  # 핫스팟 파일
panel_bed="/home/ricky8419/HMF_pipeline/resource/38/variants/ActionableCodingPanel.38.bed.gz"  # Panel BED 파일
high_confidence_bed="/home/ricky8419/HMF_pipeline/resource/38/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz"  # High Confidence BED 파일
driver_gene_panel="/home/ricky8419/HMF_pipeline/resource/38/common/DriverGenePanel.38.tsv"  # Driver Gene Panel 파일
mappability_bed="/home/ricky8419/HMF_pipeline/resource/38/variants/mappability_150.38.bed.gz"  # Mappability BED 파일
clinvar_vcf="/home/ricky8419/HMF_pipeline/resource/38/variants/clinvar.38.vcf.gz"  # ClinVar VCF 파일
blacklist_bed="/home/ricky8419/HMF_pipeline/resource/38/variants/KnownBlacklist.germline.38.bed"  # Blacklist BED 파일
blacklist_vcf="/home/ricky8419/HMF_pipeline/resource/38/variants/KnownBlacklist.germline.38.vcf.gz"  # Blacklist VCF 파일
threads=32  # 스레드 수
max_memory=256  # Java 메모리 (GB)


# .bam 파일 처리
for bam_file in "${bam_dir}"/*.bam; do
  # 샘플 ID 추출
  reference_id=$(basename "${bam_file}" .bam)
  reference_bam="${bam_file}"
  sage_output_vcf="${output_dir}/${reference_id}.sage.germline.vcf.gz"
  pave_output_vcf="${pave_output_dir}/${reference_id}.pave.germline.vcf.gz"

  # SAGE 실행
  if [[ -f "${sage_output_vcf}" ]]; then
    echo "Skipping SAGE germline for ${reference_id}, since VCF ${sage_output_vcf} exists"
  else
    echo "Running SAGE germline for ${reference_id}..."
    java -Xmx${max_memory}G -jar "${sage_jar}" \
      -tumor "${reference_id}" \
      -tumor_bam "${reference_bam}" \
      -hotspots "${hotspots}" \
      -panel_bed "${panel_bed}" \
      -high_confidence_bed "${high_confidence_bed}" \
      -ref_genome "${ref_genome}" \
      -ref_genome_version "${ref_genome_version}" \
      -ensembl_data_dir "${ensembl_dir}" \
      -ref_sample_count 0 -panel_only \
      -threads "${threads}" \
      -output_vcf "${sage_output_vcf}"

    if [[ $? -eq 0 ]]; then
      echo "SAGE germline completed successfully for ${reference_id}. Output: ${sage_output_vcf}"
    else
      echo "Error running SAGE germline for ${reference_id}"
      exit 1
    fi
  fi

  # PAVE 실행
  if [[ -f "${pave_output_vcf}" ]]; then
    echo "Skipping PAVE for ${reference_id}, since VCF ${pave_output_vcf} exists"
  else
    echo "Running PAVE for ${reference_id}..."
    java -Xmx${max_memory}G -jar "${pave_jar}" \
      -sample "${reference_id}" \
      -vcf_file "${sage_output_vcf}" \
      -ensembl_data_dir "${ensembl_dir}" \
      -driver_gene_panel "${driver_gene_panel}" \
      -ref_genome "${ref_genome}" \
      -ref_genome_version "${ref_genome_version}" \
      -mappability_bed "${mappability_bed}" \
      -clinvar_vcf "${clinvar_vcf}" \
      -blacklist_bed "${blacklist_bed}" \
      -blacklist_vcf "${blacklist_vcf}" \
      -read_pass_only \
      -output_dir "${pave_output_dir}" \
      -output_vcf_file	 "${pave_output_vcf}"
    if [[ $? -eq 0 ]]; then
      echo "PAVE completed successfully for ${reference_id}. Output directory: ${pave_output_dir}"
    else
      echo "Error running PAVE for ${reference_id}"
      exit 1
    fi
  fi
done

echo "All samples processed successfully."
