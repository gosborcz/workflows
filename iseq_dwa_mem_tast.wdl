workflow fq_bwa_mem_workflow {

  meta {
    keywords: '{"keywords": ["some", "keywords"]}'
    name: 'fq_bwa_mem'
    author: 'https://gitlab.com/MateuszMarynowski'
    copyright: 'Copyright 2019 Intelliseq'
    description: 'Generic text for task'
    changes: '{"latest": "no changes"}'

    input_fastq_1: '{"name": "(bgzipped) FASTQ 1", "type": "File", "constraints": {"extension": ["_1.fq.gz"]}, "description": "First FASTQ file (bgzipped)."}'
    input_fastq_2: '{"name": "(bgzipped) FASTQ 2", "type": "File", "constraints": {"extension": ["_2.fq.gz"]}, "description": "Second FASTQ file (bgzipped)."}'
    input_sample_id: '{"name": "Sample ID", "type": "String", "description": "Sample identifier."}'
    input_reference_genome: '{"name": "Reference genome", "type": "String", "constraints": {"values": ["hg38", "grch38-no-alt"]}, "description": "Version of the reference genome to align reads to."}'
    input_chromosome: '{"name": "Chromosome", "type": "String", "required": "false", "constraints": {"values": ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY-and-the-rest"]}, "description": "Restrict reference genome to specified chromosoms (read will be simulated only from this chromosome)."}'
    input_no_threads: '{"name": "(Advanced) Number of threads", "type": "Int", "default": "8", "description": "(Advanced) Number of threads to run BWA MEM and samtools sort."}'
    input_bwa_mem_arguments: '{"name": "(Advanced) BWA MEM arguments", "type": "String", "required": "false", "default":"-K 100000000 -v 3 ", "description": "(Advanced) Overwrite BWA MEM arguments. Arguments \'-t number of threads\' and \'-R read group header line\' cannot be overwritten by this option."}'

    output_stdout_log: '{"name": "Standard out", "type": "File", "copy": "True", "description": "Standard out"}'
    output_stderr_err: '{"name": "Standard err", "type": "File", "copy": "True", "description": "Standard error"}'
    output_bco: '{"name": "Biocompute object", "type": "File", "copy": "True", "description": "Biocompute object"}'

  }

  call fq_bwa_mem

}

task fq_bwa_mem {

  File fastq_1
  File fastq_2

  String sample_id = 'no_id_provided'

  String? chromosome
  String chromosome_defined = select_first([chromosome, ""])
  String reference_genome = "hg38" # another option: String reference_genome = "grch38-no-alt"
  String reference_genome_whole_name = if (reference_genome == "hg38") then "broad-institute-" + reference_genome else reference_genome + "-analysis-set"
  String reference_genome_scope = if defined(chromosome) then reference_genome_whole_name + "-" + chromosome_defined else reference_genome_whole_name

  Int? index

  # Tools runtime settings, paths etc.
  String reference_genome_fasta_gz = "/resources/reference-genomes/" + reference_genome_scope + "/" + reference_genome_scope + ".fa.gz"
  String reference_genome_fasta = "/resources/reference-genomes/" + reference_genome_scope + "/" + reference_genome_scope + ".fa"
  String reference_genome_fasta_fai = "/resources/reference-genomes/" + reference_genome_scope + "/" + reference_genome_scope + ".fa.fai"
  Int no_threads = 8
  String bwa_mem_arguments = "-K 100000000 -v 3 -Y "
  String RG_PL="Illumina"

  String task_name = "fq_bwa_mem"
  String task_name_with_index = if defined(index) then task_name + "_" + index else task_name
  String task_version = "latest"
  String docker_image = if defined(chromosome) then "intelliseqngs/bwa-mem-" + reference_genome + "-chr-wise:" + "2.0.0-" + chromosome_defined else "intelliseqngs/bwa-mem-" + reference_genome + ":2.0.0"

  String dollar = "$"

  command <<<
    task_name="${task_name}"; task_name_with_index="${task_name_with_index}"; task_version="${task_version}"; task_docker="${docker_image}"
    source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/after-start.sh)

    # prepare basename for bam file (intesection of basename1 and basename2)
    basename1=$( basename "${fastq_1}" ); basename2=$( basename "${fastq_2}" )
    fn=''; n=0; while [[ ${dollar}{basename1:n:1} == ${dollar}{basename2:n:1} ]]; do fn+=${dollar}{basename1:n:1}; ((n++)); done
    filename=${dollar}{fn::-1}

    # unbgzip reference genome...
    gunzip ${reference_genome_fasta_gz}

    # set readgroups
    zcat -f ${fastq_1} | head -1 | cut -d ':' -f 3,4 | sed 's/:/\./g' > rg_id
    RG_ID=`cat rg_id`
    RG_PU="$RG_ID"".""${sample_id}"
    RG_LB="${sample_id}"".library"
    RG_SM="${sample_id}"
    RG_PL="${RG_PL}"

    set -e pipefail

    # run BWA MEM, samblaster and samtools
    bwa mem \
      -t ${no_threads} \
      -R "@RG\tID:""$RG_ID""\tPU:""$RG_PU""\tPL:${RG_PL}\tLB:""$RG_LB""\tSM:""$RG_SM" \
      ${bwa_mem_arguments} \
      ${reference_genome_fasta} \
      ${fastq_1} ${fastq_2} \
        | samblaster > align-markdup.bam

    java -jar /usr/bin/picard.jar \
       SortSam \
        I=align-markdup.bam \
        O=${sample_id}-$filename"_markdup.bam" \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT

  source <(curl -s https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.0.1/before-finish.sh)
  >>>

  runtime {

    docker: docker_image
    memory: "10G"
    cpu: "1"
    maxRetries: 2

  }

  output {

    File lane_markdup_bam = glob("${sample_id}-*_markdup.bam")[0]
    File lane_markdup_bam_bai = glob("${sample_id}-*_markdup.bai")[0]

    File stdout_log = stdout()
    File stderr_log = stderr()
    File bco = "bco.json"

  }

}
