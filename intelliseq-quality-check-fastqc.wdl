# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Name: FastQC quality check
#  Authors:
#    - https://gitlab.com/dzesikahoinkis
#    - https://gitlab.com/marysiaa
#    - https://gitlab.com/marpiech
#  Copyright: Copyright 2019 Intelliseq
#  Description: >
#    Runs FastQC and generates quality check reports.
#  Changes:
#    latest:
#      - no changes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


workflow quality_check_fastqc_workflow { call quality_check_fastqc {}}

task quality_check_fastqc {

  File fastq_1
  File fastq_2

  # Inputs with defaults
  String sample_id = "no_id_provided"

  String task_name = "quality_check_fastqc"
  String task_version = "latest"
  String docker_image = "intelliseqngs/fastqc:v0.1"

  String filename = basename(fastq_1, "_1.fq.gz")

  command <<<

    ### bioobject
    if [ -e "/resources.bioobject.json" ]; then RESOURCES=$(cat /resources.bioobject.json); else RESOURCES="[]"; fi
    if [ -e "/tools.bioobject.json" ]; then RESOURCES=$(cat /tools.bioobject.json); else TOOLS="[]"; fi

    printf "{\
      \"task-name\":\"${task_name}\",\
      \"task-version\":\"${task_version}\",\
      \"docker-image\":\"${docker_image}\",\
      \"resources\":$RESOURCES,\
      \"tools\":$TOOLS\
      }" | sed 's/ //g' > bioobject.json

    fastqc ${fastq_1} -o . &> ${sample_id}.${filename}_1.fastqc.stdout.stderr.log
    fastqc ${fastq_2} -o . &> ${sample_id}.${filename}_2.fastqc.stdout.stderr.log

    unzip ${filename}_1_fastqc.zip
    unzip ${filename}_2_fastqc.zip

    cat ${filename}_1_fastqc/summary.txt | awk -F '\t' 'BEGIN {print "statistic"} {print $2}' > statistic
    cat ${filename}_1_fastqc/summary.txt | awk -F '\t' 'BEGIN {print "fastq_1-status"} {print $1}' > fastq_1-status
    cat ${filename}_2_fastqc/summary.txt | awk -F '\t' 'BEGIN {print "fastq_2-status"} {print $1}' > fastq_2-status

    paste -d '\t' statistic fastq_1-status fastq_2-status > ${sample_id}.${filename}.fastqc-statistics.txt

    cp ${filename}_1_fastqc/Images/per_base_quality.png ${sample_id}.${filename}_1.fastqc-report.first-image.png
    cp ${filename}_2_fastqc/Images/per_base_quality.png ${sample_id}.${filename}_2.fastqc-report.first-image.png

    cp ${filename}_1_fastqc/fastqc_data.txt ${sample_id}.${filename}_1.fastqc-data.txt
    cp ${filename}_2_fastqc/fastqc_data.txt ${sample_id}.${filename}_2.fastqc-data.txt

    mv ${filename}_1_fastqc.html ${sample_id}.${filename}_1.fastqc-report.html
    mv ${filename}_2_fastqc.html ${sample_id}.${filename}_2.fastqc-report.html

    mv ${filename}_1_fastqc.zip ${sample_id}.${filename}_1.fastqc-report.zip
    mv ${filename}_2_fastqc.zip ${sample_id}.${filename}_2.fastqc-report.zip

  >>>

  runtime {

    docker: docker_image
    memory: "4G"
    cpu: "1"

  }

  output {

    File fastqc_1_report_html = "${sample_id}.${filename}_1.fastqc-report.html"
    File fastqc_2_report_html = "${sample_id}.${filename}_2.fastqc-report.html"

    File fastqc_1_report_zip = "${sample_id}.${filename}_1.fastqc-report.zip"
    File fastqc_2_report_zip = "${sample_id}.${filename}_2.fastqc-report.zip"

    File fastqc_statistics_txt = "${sample_id}.${filename}.fastqc-statistics.txt"

    File fastqc_1_report_data = "${sample_id}.${filename}_1.fastqc-data.txt"
    File fastqc_2_report_data = "${sample_id}.${filename}_2.fastqc-data.txt"

    File fastqc_1_report_image = "${sample_id}.${filename}_1.fastqc-report.first-image.png"
    File fastqc_2_report_image = "${sample_id}.${filename}_2.fastqc-report.first-image.png"

    # Logs
    File fastqc_1_stdout_stderr_log = "${sample_id}.${filename}_1.fastqc.stdout.stderr.log"
    File fastqc_2_stdout_stderr_log = "${sample_id}.${filename}_2.fastqc.stdout.stderr.log"

    # @Output(required=true,directory="/quality_check_fastqc",filename="stdout.log")
    File stdout_log = stdout()
    # @Output(required=true,directory="/quality_check_fastqc",filename="stderr.log")
    File stderr_log = stderr()
    # @Output(required=true,directory="/quality_check_fastqc",filename="bioobject.json")
    File bioobject = "bioobject.json"

  }

}
