# Germline Targeted Sequencing Pipeline

import "https://gitlab.com/intelliseq/workflows/raw/fq-organize@1.0.0/src/main/wdl/tasks/fq-organize/fq-organize.wdl" as fq_organize_task
import "https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/resources-kit/latest/resources-kit.wdl" as resources_kit_task
import "https://gitlab.com/intelliseq/workflows/raw/fq-qc@1.0.1/src/main/wdl/modules/fq-qc/latest/fq-qc.wdl" as fq_qc_module
import "https://gitlab.com/intelliseq/workflows/raw/alignment@1.2.3/src/main/wdl/modules/alignment/latest/alignment.wdl" as alignment_module
import "https://gitlab.com/intelliseq/workflows/raw/variant-calling@1.0.2/src/main/wdl/modules/variant-calling/latest/variant-calling.wdl" as variant_calling_module
import "https://gitlab.com/intelliseq/workflows/raw/vcf-anno@1.0.0/src/main/wdl/modules/vcf-anno/latest/vcf-anno.wdl" as vcf_anno_module
import "https://gitlab.com/intelliseq/workflows/raw/vcf-acmg-report@1.0.3/src/main/wdl/modules/vcf-acmg-report/latest/vcf-acmg-report.wdl" as vcf_acmg_report_module
import "https://gitlab.com/intelliseq/workflows/raw/coverage-statistics@1.0.0/src/main/wdl/modules/coverage-statistics/latest/coverage-statistics.wdl" as coverage_statistics_module
import "https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/merge-pdf/latest/merge-pdf.wdl" as merge_pdf_task
import "https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/merge-bco-pipeline/latest/merge-bco-pipeline.wdl" as merge_bco_pipeline_task
import "https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/report-bco/latest/report-bco.wdl" as report_bco_task
#import "https://gitlab.com/intelliseq/workflows/raw/detection-chance@1.1.1/src/main/wdl/modules/detection-chance/latest/detection-chance.wdl" as detection_chance_module


workflow target_seq {


  meta {
    name: 'Germline target-seq workflow'
    author: 'https://gitlab.com/marpiech'
    copyright: 'Copyright 2019-2020 Intelliseq'
    description: '## Workflow for identification of pathogenic variants\nValidated for gene panels (including hereditary cancer) and whole exome data.  Generates diagnostic report. CE certified pending.'
    released: 'true'
    input_fastqs: '{"name": "fastq files", "type": "Array[File]", "constraints": {"extension": ["fq.gz"]}, "description": "List of fastq files (both left and right)"}'
    input_fastqs_left: '{"hidden":"true", "name": "fastq 1", "type": "Array[File]", "constraints": {"extension": ["fq.gz"]}, "description": "first fastq file"}'
    input_fastqs_right: '{"hidden":"true", "name": "fastq 2", "type": "Array[File]", "constraints": {"extension": ["fq.gz"]}, "description": "second fastq file"}'
    input_sample_id: '{"name": "sample id", "type": "String", "description": "identifier of sample"}'
    input_genome_or_exome: '{"name": "genome_or_exome", "type": "String", "default": "exome", "description": "Gnomad coverage"}'
    input_interval_list: '{"name": "Intervals",  "type": "File", "constraints": {"extension": "interval_list"}, "description": "List of genomic intervals. Default sure-select exome v7"}'
    input_panel_json: '{"name": "panel_json", "type": "File",  "constraints": {"extension": ["json"]}, "description": "Panel data in json format, must include attributes (even if without values): genes_number, genes"}'
    input_patient_json: '{"name": "patient_json", "type": "File",  "constraints": {"extension": ["json"]}, "description": "Patient data in json format, must include attributes (even if without values): name, surname, sex, birthdate and pesel"}'
    input_sample_json: '{"name": "sample_json", "type": "File",  "constraints": {"extension": ["json"]}, "description": "Patient data in json format, must include attributes (even if without values): ID, material, sequencing_type, sequencing_platform, sending_date, raport_date, doctor_name "}'
    input_phenotypes_json: '{"name": "phenotypes_json", "type": "File",  "constraints": {"extension": ["json"]}, "description": "Contains info about used phenotypes with HPid number, list of genes and number of genes used to call"}'
    input_other_bams: '{"name": "Other bams", "type": "Array[File]",  "constraints": {"extension": ["bam"]}, "description": "Bams files used for igv-screenshots comparison"}'
    input_other_bais: '{"name": "Other bais", "type": "Array[File]",  "constraints": {"extension": ["bai"]}, "description": "Bais files needed for igv-screenshots comparison"}'

    output_final_bam: '{"name": "bam", "type": "File", "copy": "True", "description": "Alignment result"}'
    output_final_bai: '{"name": "bai", "type": "File", "copy": "True", "description": "Alignment result index"}'
    output_gvcf_gz: '{"name": "g.vcf.gz", "type": "File", "copy": "True", "description": "Genotyping result"}'
    output_gvcf_gz_tbi: '{"name": "g.vcf.gz.tbi", "type": "File", "copy": "True", "description": "Genotyping result index"}'
    output_vcf_gz: '{"name": "vcf.gz", "type": "File", "copy": "True", "description": "Variants"}'
    output_vcf_gz_tbi: '{"name": "vcf.gz.tbi", "type": "File", "copy": "True", "description": "Variants index"}'
    output_annotated_vcf: '{"name": "annotated.vcf.gz", "type": "File", "copy": "True", "description": "Annotation result"}'
    output_annotated_vcf_tbi: '{"name": "annotated.vcf.gz.tbi", "type": "File", "copy": "True", "description": "Annotation result index"}'
    output_ang_pdf_report: '{"name": "Report from genetic analysis in English", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, pdf format, English version"}'
    output_ang_html_report: '{"name": "Report from genetic analysis in English", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, html format, English version"}'
    output_ang_docx_report: '{"name": "Report from genetic analysis in English", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, docx format, English version"}'
    output_ang_odt_report: '{"name": "Report from genetic analysis in English", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, odt format, English version"}'
    output_pdf_report: '{"name": "Report from genetic analysis in Polish", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, pdf format, Polish version"}'
    output_html_report: '{"name": "Report from genetic analysis in Polish", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, html format, Polish version"}'
    output_docx_report: '{"name": "Report from genetic analysis in Polish", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, docx format, Polish version"}'
    output_odt_report: '{"name": "Report from genetic analysis in Polish", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, odt format, Polish version"}'
    output_csv: '{"name": "csv", "type": "File", "copy": "True", "description": "Table with variants"}'
    output_fastqc_1_zips: '{"name": "FastQC output zip", "type": "Array[File]", "copy": "True", "description": "FastQC results for the first fastq files"}'
    output_fastqc_2_zips: '{"name": "FastQC output zip", "type": "Array[File]", "copy": "True", "description": "FastQC results for the second fastq files"}'
    output_quality_check_1_jsons: '{"name": "Fastq quality check json", "type": "Array[File]", "copy": "True", "description": "Jsons with quality check results for the first fastq files"}'
    output_quality_check_2_jsons: '{"name": "Fastq quality check json", "type": "Array[File]", "copy": "True", "description": "Jsons with quality check results for the second fastq files"}'
    output_quality_check_pdf: '{"name": "Quality check pdf", "type": "File", "copy": "True", "description": "Fastq quality check report in pdf format"}'
    output_simple_coverage_stats_report_pdf: '{"name": "pdf", "type": "File", "copy": "True", "description": "pdf human readeble coverage statistics report"}'
    output_coverage_stats_report_pdf: '{"name": "pdf", "type": "File", "copy": "True", "description": "pdf coverage statistics report based on gatk nad samtools"}'
    output_full_report_pdf: '{"name": "pdf", "type": "File", "copy": "True", "description": "all pdf reports merged into one pdf file"}'
    output_final_realigned_bam: '{"name": "bam", "type": "File", "copy": "True", "description": "Alignment result"}'
    output_final_realigned_bai: '{"name": "bai", "type": "File", "copy": "True", "description": "Alignment result index"}'
    output_igv_screenshots_tar_gz: '{"name": "igv screenshots tar", "type": "File", "copy": "True", "description": "Igv screenshots of alignment at variant positions"}'
  }


  String sample_id = "no_id_provided"
  Array[File]? fastqs
  Boolean is_fastqs_defined = defined(fastqs)
  Array[File]? fastqs_left
  Array[File]? fastqs_right
  File? panel_json
  File? patient_json
  File? sample_json
  File? phenotypes_json
  File? interval_list
  File? error_warning_json_qc
  File? error_warning_json_coverage
  String genome_or_exome = "exome"
  String kit = "exome-v7"
  Boolean is_interval_list_not_defined = !defined(interval_list)
  Array[File]? other_bams = []
  Array[File]? other_bais = []
  String pipeline_name = "target_seq"
  String pipeline_version = "latest"

  ### organizes fastqs
  if(is_fastqs_defined) {
    call fq_organize_task.fq_organize {
      input:
      fastqs = fastqs
    }
  }
  Array[File] fastqs_1 = select_first([fq_organize.fastqs_1, fastqs_left])
  Array[File] fastqs_2 = select_first([fq_organize.fastqs_2, fastqs_left])

  ### fills interval_list if not declared ###
  if(is_interval_list_not_defined) {
    call resources_kit_task.resources_kit {
      input:
      kit = kit
    }
  }
  File interval_file = select_first([interval_list, resources_kit.interval_list])

  ### checks quality of fastq files ###
  call fq_qc_module.fq_qc {
    input:
      sample_id = sample_id,
      fastqs_1 = fastqs_1,
      fastqs_2 = fastqs_2
  }

  call alignment_module.alignment {
    input:
      fastqs_1 = fastqs_1,
      fastqs_2 = fastqs_2,
      sample_id = sample_id
  }

  call variant_calling_module.variant_calling {
    input:
      input_bam = alignment.recalibrated_markdup_bam,
      input_bai = alignment.recalibrated_markdup_bai,
      sample_id = sample_id,
      interval_list = interval_file
  }

  call vcf_anno_module.vcf_anno {
    input:
      vcf_gz = variant_calling.vcf_gz,
      vcf_gz_tbi = variant_calling.vcf_gz_tbi,
      gnomad_coverage_genome_or_exome = genome_or_exome,
      vcf_basename = sample_id
  }

  #call detection_chance_module.detection_chance {
  #  input:
  #    sample_gvcf_gz = variant_calling.gvcf_gz,
  #    sample_gvcf_gz_tbi = variant_calling.gvcf_gz_tbi,
  #    panel_json = panel_json
  #}

  call vcf_acmg_report_module.vcf_acmg_report {
    input:
      panel_json = panel_json,
      sample_json = sample_json,
      phenotypes_json = phenotypes_json,
      patient_json = patient_json,
      vcf_gz = vcf_anno.annotated_and_filtered_vcf,
      vcf_gz_tbi = vcf_anno.annotated_and_filtered_vcf_tbi,
      bam = alignment.recalibrated_markdup_bam,
      bai = alignment.recalibrated_markdup_bai,
      realigned_bam = variant_calling.haplotype_caller_bam,
      realigned_bai = variant_calling.haplotype_caller_bai,
      other_bams = other_bams,
      other_bais = other_bais,
      sample_id = sample_id
  }

#  call coverage_statistics_module.coverage_statistics {
#    input:
#      intervals = interval_file,
#      bam_file = alignment.recalibrated_markdup_bam,
#      bai_file = alignment.recalibrated_markdup_bai,
#      error_warning_json_coverage = error_warning_json_coverage
#  }

#  Array[File] reports_pdf = [vcf_acmg_report.pdf_report,fq_qc.quality_check_report_attach_pdf,coverage_statistics.simple_coverage_report_pdf]
  Array[File] reports_pdf = [vcf_acmg_report.pdf_report,fq_qc.quality_check_report_pdf]
#  Array[File] reports_ang_pdf = [vcf_acmg_report.ang_pdf_report,fq_qc.quality_check_report_attach_pdf,coverage_statistics.simple_coverage_report_pdf]
  Array[File] reports_ang_pdf = [vcf_acmg_report.ang_pdf_report,fq_qc.quality_check_report_pdf]

  call merge_pdf_task.merge_pdf {
    input:
      sample_id = sample_id,
      pdf = reports_pdf,
      ang_pdf = reports_ang_pdf
  }

  Array[File] bcos_module = [fq_qc.bco, alignment.bco, variant_calling.bco, vcf_anno.bco]
  Array[File] stdout_module = [fq_qc.stdout_log, alignment.stdout_log, variant_calling.stdout_log, vcf_anno.stdout_log]
  Array[File] stderr_module = [fq_qc.stderr_log, alignment.stderr_log, variant_calling.stderr_log, vcf_anno.stderr_log]

  call merge_bco_pipeline_task.merge_bco_pipeline {
    input:
        bcos = bcos_module,
        stdout = stdout_module,
        stderr = stderr_module,
        pipeline_name = pipeline_name,
        pipeline_version = pipeline_version
  }

  call report_bco_task.report_bco as report_bco_pipeline {
    input:
       bco_json = merge_bco_pipeline.bco
  }

  output {
    File final_bam = alignment.recalibrated_markdup_bam
    File final_bai = alignment.recalibrated_markdup_bai

    File final_realigned_bam = variant_calling.haplotype_caller_bam
    File final_realigned_bai = variant_calling.haplotype_caller_bai

    File gvcf_gz = variant_calling.gvcf_gz
    File gvcf_gz_tbi = variant_calling.gvcf_gz_tbi

    File vcf_gz = variant_calling.vcf_gz
    File vcf_gz_tbi = variant_calling.vcf_gz_tbi

    File annotated_vcf = vcf_anno.annotated_and_filtered_vcf
    File annotated_vcf_tbi = vcf_anno.annotated_and_filtered_vcf_tbi

    #report
    File csv_report = vcf_acmg_report.csv_report
    File annotated_acmg_vcf_gz = vcf_acmg_report.annotated_acmg_vcf_gz
    File annotated_acmg_vcf_gz_tbi = vcf_acmg_report.annotated_acmg_vcf_gz_tbi
    File docx_report = vcf_acmg_report.docx_report
    File pdf_report = vcf_acmg_report.pdf_report
    File html_report = vcf_acmg_report.html_report
    File odt_report = vcf_acmg_report.odt_report
    File ang_docx_repor = vcf_acmg_report.ang_docx_repor
    File ang_pdf_report = vcf_acmg_report.ang_pdf_report
    File ang_html_report = vcf_acmg_report.ang_html_report
    File ang_odt_report = vcf_acmg_report.ang_odt_report


    File bco_report_pdf = report_bco_pipeline.bco_report_pdf
    File bco_report_odt = report_bco_pipeline.bco_report_odt
    File bco_report_docx = report_bco_pipeline.bco_report_docx
    File bco_report_html = report_bco_pipeline.bco_report_html

    #File detection_chance_report_pdf = detection_chance.detection_chance_report_pdf
    #File detection_chance_report_odt = detection_chance.detection_chance_report_odt
    #File detection_chance_report_docx = detection_chance.detection_chance_report_docx
    #File detection_chance_report_html = detection_chance.detection_chance_report_html

    #qc
    Array[File] fastqc_1_zips = fq_qc.fastqc_1_zips
    Array[File] fastqc_2_zips = fq_qc.fastqc_2_zips
    Array[File] quality_check_1_jsons = fq_qc.quality_check_1_jsons
    Array[File] quality_check_2_jsons = fq_qc.quality_check_2_jsons
    File quality_check_pdf = fq_qc.quality_check_report_pdf

#    #coverage statistics
#    File coverage_stats_report_pdf = coverage_statistics.coverage_report_pdf
#    File simple_coverage_stats_report_pdf = coverage_statistics.simple_coverage_report_pdf

    #full report
    File full_report_pdf = merge_pdf.full_report_pdf
    File full_report_ang_pdf = merge_pdf.full_report_ang_pdf

    #igv-screenshots
    File igv_screenshots_tar_gz = vcf_acmg_report.igv_screenshots_tar_gz

    #bco, stdout
    File bco = merge_bco_pipeline.bco
    File stdout_log = merge_bco_pipeline.stdout_log

  }
}
