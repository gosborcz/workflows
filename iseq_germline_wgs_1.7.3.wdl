import "https://gitlab.com/intelliseq/workflows/raw/fq-organize@1.1.4/src/main/wdl/tasks/fq-organize/fq-organize.wdl" as fq_organize_task
import "https://gitlab.com/intelliseq/workflows/raw/panel-generate@1.6.5/src/main/wdl/tasks/panel-generate/panel-generate.wdl" as panel_generate_task
import "https://gitlab.com/intelliseq/workflows/raw/resources-kit@1.0.1/src/main/wdl/tasks/resources-kit/resources-kit.wdl" as resources_kit_task
import "https://gitlab.com/intelliseq/workflows/raw/fq-qc@1.4.1/src/main/wdl/modules/fq-qc/fq-qc.wdl" as fq_qc_module
import "https://gitlab.com/intelliseq/workflows/raw/fq-bwa-align@1.4.1/src/main/wdl/modules/fq-bwa-align/latest/fq-bwa-align.wdl" as fq_bwa_align_module
import "https://gitlab.com/intelliseq/workflows/raw/bam-filter-contam@1.0.4/src/main/wdl/modules/bam-filter-contam/bam-filter-contam.wdl" as bam_filter_contam_module
import "https://gitlab.com/intelliseq/workflows/raw/sv-calling@1.0.8/src/main/wdl/modules/sv-calling/sv-calling.wdl" as sv_calling_module
import "https://gitlab.com/intelliseq/workflows/raw/bam-varcalling@1.3.0/src/main/wdl/modules/bam-varcalling/bam-varcalling.wdl" as bam_varcalling_module
import "https://gitlab.com/intelliseq/workflows/raw/vcf-anno@1.7.1/src/main/wdl/modules/vcf-anno/vcf-anno.wdl" as vcf_anno_module
import "https://gitlab.com/intelliseq/workflows/raw/vcf-acmg-report@1.1.9/src/main/wdl/modules/vcf-acmg-report/latest/vcf-acmg-report.wdl" as vcf_acmg_report_module
import "https://gitlab.com/intelliseq/workflows/raw/bam-qc@2.2.2/src/main/wdl/modules/bam-qc/bam-qc.wdl" as bam_qc_module
import "https://gitlab.com/intelliseq/workflows/raw/detection-chance@1.3.9/src/main/wdl/modules/detection-chance/latest/detection-chance.wdl" as detection_chance_module
import "https://gitlab.com/intelliseq/workflows/raw/sex-check@1.0.2/src/main/wdl/modules/sex-check/sex-check.wdl" as sex_check_module
import "https://gitlab.com/intelliseq/workflows/raw/vcf-var-filter@1.0.2/src/main/wdl/modules/vcf-var-filter/vcf-var-filter.wdl" as vcf_var_filter_module
import "https://gitlab.com/intelliseq/workflows/raw/pdf-merge@1.1.1/src/main/wdl/tasks/pdf-merge/latest/pdf-merge.wdl" as pdf_merge_task
import "https://gitlab.com/intelliseq/workflows/raw/bco-merge@1.4.0/src/main/wdl/tasks/bco-merge/latest/bco-merge.wdl" as bco_merge_task
import "https://gitlab.com/intelliseq/workflows/raw/report-bco@1.0.1/src/main/wdl/tasks/report-bco/latest/report-bco.wdl" as report_bco_task

workflow germline {

    meta {
        author: 'https://gitlab.com/marpiech'
        copyright: 'Copyright 2019-2020 Intelliseq'
        tag: 'Clinical WES/WGS'

        groups: '{"gene_panel": {"description": "Fill in at least 1 input below to generate the gene panel", "min_inputs": 1}}'

        variant1_name: 'WES hereditary disorders ACMG report'
        variant1_description: 'Identifies pathogenic variants (according to the ACMG classification) with the ready-to-use or custom gene panels; generates full clinical report; for whole exome sequencing data (WES).'

        variant1_input_genes: '{"index": 5, "name": "Genes names", "type": "String", "groupname": "gene_panel", "description": "Enter gene names to narrow your search/analysis results (separate gene names with comma, for example: HTT, FBN1)"}'
        variant1_input_genome_or_exome: '{"hidden":"true", "name": "Sample type", "type": "String", "default": "exome", "constraints": {"values": ["exome", "genome"]}, "description": "Select exome or genome (refers to the input sample)"}'
        variant1_input_run_sv_calling: '{"hidden":"true", "name": "Run SV calling", "type": "Boolean", "default": false, "description": "Select if run SV calling"}'
        variant1_input_kit: '{"index": 8, "advanced":"true", "name": "Kit", "type": "Array[String]", "default": "exome-v7", "constraints": {"values": ["exome-v6", "exome-v7"], "multiselect": false}, "description": "Choose comprehensive exome kit from exome-v6 for Agilent SureSelect Human All Exon V6 or exome-v7 for Agilent SureSelect Human All Exon V7"}'
        variant1_input_interval_list: '{"index": 9, "hidden":"true", "name": "Intervals", "type": "File", "extension": [".interval_list"], "description": "List of genomic intervals"}'

        variant2_name: 'WES demo ACMG report'
        variant2_description: 'Demo: identifies pathogenic variants (according to the ACMG classification); generates full clinical report; for whole exome sequencing data (WES).'

        variant2_input_genes: '{"index": 5, "name": "Genes names", "type": "String", "groupname": "gene_panel", "default": "GNE", "description": "Enter gene names to narrow your search/analysis results (separate gene names with comma, for example: HTT, FBN1)"}'
        variant2_input_genome_or_exome: '{"hidden":"true", "name": "Sample type", "type": "String", "default": "exome", "constraints": {"values": ["exome", "genome"]}, "description": "Select exome or genome (refers to the input sample)"}'
        variant2_input_run_sv_calling: '{"hidden":"true", "name": "Run SV calling", "type": "Boolean", "default": false, "description": "Select if run SV callung"}'
        variant2_input_kit: '{"index": 8, "advanced":"true", "name": "Kit", "type": "Array[String]", "default": "exome-v7", "constraints": {"values": ["exome-v6", "exome-v7"], "multiselect": false}, "description": "Choose comprehensive exome kit from exome-v6 for Agilent SureSelect Human All Exon V6 or exome-v7 for Agilent SureSelect Human All Exon V7"}'
        variant2_input_interval_list: '{"index": 9, "hidden":"true", "name": "Intervals", "type": "File", "extension": ["interval_list"], "description": "List of genomic intervals"}'

        variant3_name: 'WGS hereditary disorders ACMG report'
        variant3_description: 'Identifies pathogenic variants (according to the ACMG classification) with the ready-to-use or custom gene panels; generates full clinical report; for whole genome sequencing data (WGS).'

        variant3_input_genes: '{"index": 5, "name": "Genes names", "type": "String", "groupname": "gene_panel", "description": "Enter gene names to narrow your search/analysis results (separate gene names with comma, for example: HTT, FBN1)"}'
        variant3_input_genome_or_exome: '{"hidden":"true", "name": "Sample type", "type": "String", "default": "genome", "constraints": {"values": ["exome", "genome"]}, "description": "Select exome or genome (refers to the input sample)"}'
        variant3_input_run_sv_calling: '{"hidden":"true", "name": "Run SV calling", "type": "Boolean", "default": false, "description": "Select if run SV calling"}'
        variant3_input_kit: '{"index": 8, "hidden":"true", "name": "Kit", "type": "Array[String]", "default": "genome", "constraints": {"values": ["genome", "exome-v6", "exome-v7"], "multiselect": false}, "description": "Choose comprehensive exome kit from exome-v6 for Agilent SureSelect Human All Exon V6 or exome-v7 for Agilent SureSelect Human All Exon V7"}'
        variant3_input_interval_list: '{"index": 9, "advanced":"true", "name": "Intervals", "type": "File", "extension": [".interval_list"], "description": "List of genomic intervals"}'

        variant4_name: 'WGS demo ACMG report'
        variant4_description: 'Demo: identifies pathogenic variants (according to the ACMG classification); generates full clinical report; for whole genome sequencing data (WGS).'

        variant4_input_genes: '{"index": 5, "name": "Genes names", "type": "String", "groupname": "gene_panel", "default": "GNE", "description": "Enter gene names to narrow your search/analysis results (separate gene names with comma, for example: HTT, FBN1)"}'
        variant4_input_genome_or_exome: '{"hidden":"true", "name": "Sample type", "type": "String", "default": "genome", "constraints": {"values": ["exome", "genome"]}, "description": "Select exome or genome (refers to the input sample)"}'
        variant4_input_run_sv_calling: '{"hidden":"true", "name": "Run SV calling", "type": "Boolean", "default": false, "description": "Select if run SV calling"}'
        variant4_input_kit: '{"index": 8, "hidden":"true", "name": "Kit", "type": "Array[String]", "default": "genome", "constraints": {"values": ["genome", "exome-v6", "exome-v7"], "multiselect": false}, "description": "Choose comprehensive exome kit from exome-v6 for Agilent SureSelect Human All Exon V6 or exome-v7 for Agilent SureSelect Human All Exon V7"}'
        variant4_input_interval_list: '{"index": 9, "advanced":"true", "name": "Intervals", "type": "File", "extension": [".interval_list"], "description": "List of genomic intervals"}'

        input_sample_id: '{"index": 1, "name": "Sample id", "type": "String", "default": "no_id_provided", "description": "Enter a sample name (or identifier)"}'
        input_fastqs: '{"index": 2, "name": "Fastq files", "required": "true", "paired": "true", "type": "Array[File]", "extension": [".fq.gz", ".fastq.gz"], "description": "Choose list of paired gzipped fastq files both left and right [.fq.gz or .fastq.gz]"}'
        input_fastqs_left: '{"hidden":"true", "name": "First (left) fastq files", "type": "Array[File]", "extension": [".fq.gz", ".fastq.gz"], "description": "Choose first (left) fastq files"}'
        input_fastqs_right: '{"hidden":"true", "name": "Second (right) fastq files", "type": "Array[File]", "extension": [".fq.gz", ".fastq.gz"], "description": "Choose second (right) fastq files"}'

        input_hpo_terms: '{"index": 3, "name": "HPO terms", "type": "String", "groupname": "gene_panel", "description": "Enter HPO terms to narrow your search/analysis results (separate HPO terms with comma, for example: HP:0004942, HP:0011675)"}'
        input_diseases: '{"index": 4, "name": "Diseases", "type": "String", "groupname": "gene_panel","description": "Enter disease names to narrow your search/analysis results (separate diseases names with comma; each disease name should be just a keyword, for example for Marfan Syndrome only Marfan should be written, for Ehlers-Danlos Syndrome: Ehlers-Danlos; other proper diseases names for example: Osteogenesis imperfecta, Tay-sachs, Hemochromatosis, Brugada, Canavan, etc.)"}'
        input_phenotypes_description: '{"index": 12, "name": "Description of patient phenotypes", "type": "String", "groupname": "gene_panel", "description": "Enter description of patient phenotypes"}'
        input_panel_names: '{"index": 13, "name": "Gene panel", "type": "Array[String]", "groupname": "gene_panel", "description": "Select gene panels", "constraints": {"values": ["None", "ACMG_Incidental_Findings", "COVID-19_research", "Cancer_Germline", "Cardiovascular_disorders", "Ciliopathies", "Dermatological_disorders", "Dysmorphic_and_congenital_abnormality_syndromes", "Endocrine_disorders", "Gastroenterological_disorders", "Growth_disorders", "Haematological_and_immunological_disorders", "Haematological_disorders", "Hearing_and_ear_disorders", "Metabolic_disorders", "Neurology_and_neurodevelopmental_disorders", "Ophthalmological_disorders", "Rare_Diseases", "Renal_and_urinary_tract_disorders", "Respiratory_disorders", "Rheumatological_disorders", "Skeletal_disorders", "Tumour_syndromes"], "multiselect": true}}'

        input_panel_json: '{"hidden":"true", "name": "Genes panel", "type": "File", "extension": [".json"], "description": "Add json file with genes panel to narrow your search/analysis results. You can prepare this file using workflow Gene panel generator"}'
        input_phenotypes_json: '{"hidden":"true", "name": "List of phenotypes names", "type": "File", "extension": [".json"], "description": "Add json file with phenotypes names used to generate gene panel. You can prepare this file using workflow Gene panel generator"}'

        input_patient_json: '{"index": 6, "advanced":"true", "name": "Patient information", "type": "File", "extension": [".json"], "description": "Add patient data in json format, must include attributes (even if without values): name, surname, sex, birthdate and pesel"}'
        input_sample_json: '{"index": 7, "advanced":"true", "name": "Sample information ", "type": "File", "extension": [".json"], "description": "Add sample data in json format, must include attributes (even if without values): ID, material, sequencing_type, sequencing_platform, sending_date, raport_date, doctor_name"}'

        input_other_bams: '{"index": 10, "advanced":"true", "name": "Other bams", "type": "Array[File]", "extension": [".bam"], "description": "Bams files used for igv-screenshots comparison"}'
        input_other_bais: '{"index": 11, "advanced":"true", "name": "Other bais", "type": "Array[File]", "extension": [".bai"], "description": "Bais files needed for igv-screenshots comparison"}'

        input_add_bam_filter: '{"hidden":"true", "name": "Add bam filter?", "type": "Boolean", "default": "true", "description": "This option decides whether to apply bam filtering on contamination. Default: true"}'
        input_mismatch_threshold: '{"hidden":"true", "name": "Mismatch threshold", "type": "Float", "constraints": {"min": "0.0", "max": "1.0"}, "default": "0.1", "description": "Reads with fraction of mismatches greater than or equal this value are discarded"}'


        output_fastqc_1_zips: '{"name": "FastQC output zip", "type": "Array[File]", "copy": "True", "description": "FastQC results for the first fastq files"}'
        output_fastqc_2_zips: '{"name": "FastQC output zip", "type": "Array[File]", "copy": "True", "description": "FastQC results for the second fastq files"}'
        output_quality_check_pdf: '{"name": "Quality check pdf", "type": "File", "copy": "True", "description": "Fastq quality check report in pdf format"}'
        output_quality_check_html: '{"name": "Quality check html", "type": "File", "copy": "True", "description": "Fastq quality check report in html format"}'

        output_final_bam: '{"name": "bam", "type": "File", "copy": "True", "description": "Alignment result"}'
        output_final_bai: '{"name": "bai", "type": "File", "copy": "True", "description": "Alignment result index"}'

        output_gvcf_gz: '{"name": "g.vcf.gz", "type": "File", "copy": "True", "description": "Genotyping result"}'
        output_gvcf_gz_tbi: '{"name": "g.vcf.gz.tbi", "type": "File", "copy": "True", "description": "Genotyping result index"}'
        output_vcf_gz: '{"name": "vcf.gz", "type": "File", "copy": "True", "description": "Variants"}'
        output_vcf_gz_tbi: '{"name": "vcf.gz.tbi", "type": "File", "copy": "True", "description": "Variants index"}'
        output_final_realigned_bam: '{"name": "realigned_bam", "type": "File", "copy": "True", "description": "realigned_bam"}'
        output_final_realigned_bai: '{"name": "realigned_bai", "type": "File", "copy": "True", "description": "realigned_bai"}'

        output_filtered_vcf_gz: '{"name": "Filtered vcf gz", "type": "File", "copy": "True", "description": "VCF file with variants that passed all applied filters (bgzipped)"}'
        output_filtered_vcf_gz_tbi: '{"name": "Filtered vcf gz tbi", "type": "File", "copy": "True", "description": "Index for the filtered vcf file"}'
        output_rejected_variants_vcfs: '{"name": "Rejected variants vcfs gz", "type": "File[Array]", "copy": "True", "description": "VCF files with variants that failed any of filters applied (bgzipped)"}'
        output_rejected_variants_tbis: '{"name": "Rejected variants tbis", "type": "File[Array]", "copy": "True", "description": "Indexes for the rejected variants vcfs"}'
        output_detail_metrics_file: '{"name": "Detail metrics file", "type": "File", "copy": "True", "description": "File with variant calling detailed metrics"}'
        output_summary_metrics_file: '{"name": "Summary metrics file", "type": "File", "copy": "True", "description": "File with variant calling summary metrics"}'

        output_annotated_vcf: '{"name": "annotated.vcf.gz", "type": "File", "copy": "True", "description": "Annotation result"}'
        output_annotated_vcf_tbi: '{"name": "annotated.vcf.gz.tbi", "type": "File", "copy": "True", "description": "Annotation result index"}'

        output_ang_pdf_report: '{"name": "Report from genetic analysis in English", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, pdf format, English version"}'
        output_ang_html_report: '{"name": "Report from genetic analysis in English", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, html format, English version"}'
        output_ang_docx_report: '{"name": "Report from genetic analysis in English", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, docx format, English version"}'
        output_ang_odt_report: '{"name": "Report from genetic analysis in English", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, odt format, English version"}'

        output_pdf_report: '{"name": "Report from genetic analysis in Polish", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, pdf format, Polish version"}'
        output_html_report: '{"name": "Report from genetic analysis in Polish", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, html format, Polish version"}'
        output_docx_report: '{"name": "Report from genetic analysis in Polish", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, docx format, Polish version"}'
        output_all_reports_pdf: '{"name": "All pdf reports from genetic analysis", "type": "Array[File]", "copy": "True", "description": "All pdf reports from genetic analysis"}'
        output_odt_report: '{"name": "Report from genetic analysis in Polish", "type": "File", "copy": "True", "description": "Report with results of the genetic analysis, odt format, Polish version"}'
        output_csv: '{"name": "csv", "type": "File", "copy": "True", "description": "Table with variants"}'
        output_igv_screenshots_tar_gz: '{"name": "igv screenshots tar", "type": "File", "copy": "True", "description": "Igv screenshots of alignment at variant positions"}'

        output_simple_coverage_stats_report_pdf: '{"name": "pdf", "type": "File", "copy": "True", "description": "pdf human readable coverage statistics report"}'
        output_coverage_stats_report_pdf: '{"name": "pdf", "type": "File", "copy": "True", "description": "pdf coverage statistics report based on gatk nad samtools"}'
        output_full_report_pdf: '{"name": "pdf", "type": "File", "copy": "True", "description": "all pdf reports merged into one pdf file"}'
    }

    String sample_id = "no_id_provided"

    # Provide fastq files if you would like to start analysis from the begining
    Array[File]? fastqs
    Boolean is_fastqs_defined = defined(fastqs)
    Array[File]? fastqs_left
    Boolean is_fastqs_left_defined = defined(fastqs_left)
    Array[File]? fastqs_right

    # Provide bam file if you would like to start analysis from variant calling
    File? bam
    File? bam_bai
    Boolean is_bam_defined = defined(bam)

    # Gene panel
    File? panel_json
    File? phenotypes_json
    Boolean is_panel_json_defined = defined(panel_json)
    String? hpo_terms
    String? genes
    String? diseases
    String? phenotypes_description
    Array[String]? panel_names
    Boolean is_input_for_panel_generate_defined = (defined(hpo_terms) || defined(genes) || defined(diseases) || defined(phenotypes_description) || defined(panel_names))

    # Patient information
    File? patient_json
    File? sample_json

    File? error_warning_json_qc # It needs to be add to meta section
    File? error_warning_json_coverage # It needs to be add to meta section

    File? interval_list
    String genome_or_exome = "exome" # Possible values "exome", "genome"
    String kit = "exome-v7" # Possible values "exome-v6", "exome-v7", "genome"
    String kit_choice = if (genome_or_exome == "genome") then "genome" else kit
    Boolean is_interval_list_not_defined = !defined(interval_list)

    Int max_no_pieces_to_scatter_an_interval_file = 6
    Array[File]? other_bams = []
    Array[File]? other_bais = []
    Boolean add_bam_filter = true
    Float mismatch_threshold = 0.1

    # sv parameters
    Boolean run_sv_calling = false
    Boolean add_snp_data_to_sv = true
    Float del_cov_max = 0.75
    Float dup_cov_min = 1.3
    Float del_cov_gcbin_max = 0.75
    Float dup_cov_gcbin_min = 1.3
    Boolean exclude_lcr = true

    # AD filtering
    Float ad_binom_threshold = 0.01
    Boolean apply_ad_filter = false

    # IGV (make something much larger if pictures in report are not fully loaded)
    Int waiting_time = 10000

    String pipeline_name = "germline"
    String pipeline_version = "1.7.3"

    # 1. Prepare gene panel or use user defined
    if(is_input_for_panel_generate_defined && !is_panel_json_defined) {
        call panel_generate_task.panel_generate {
            input:
                sample_id = sample_id,
                hpo_terms = hpo_terms,
                genes = genes,
                diseases = diseases,
                phenotypes_description = phenotypes_description,
                panel_names = panel_names
        }
    }

    # Use acmg gene panel if panel and panels inputs are not defined
    if(!is_input_for_panel_generate_defined && !is_panel_json_defined) {
        call panel_generate_task.panel_generate as panel_generate_acmg_panel {
            input:
                sample_id = sample_id,
                panel_names = ["acmg-recommendation-panel"]
        }
    }

    File gene_panel_json = select_first([panel_json, panel_generate.panel, panel_generate_acmg_panel.panel])
    File gene_phenotypes_json = select_first([phenotypes_json, panel_generate.phenotypes, panel_generate_acmg_panel.phenotypes])

    # 2. Prepare interval_list
    if(is_interval_list_not_defined) {
        call resources_kit_task.resources_kit {
            input:
                kit = kit_choice
        }
    }

    File interval_file = select_first([interval_list, resources_kit.interval_list])

    # Start analysis form fastq file
    if(is_fastqs_defined || is_fastqs_left_defined) {
        if(is_fastqs_defined) {

            # 3. Organise fastq files
            call fq_organize_task.fq_organize {
                input:
                    fastqs = fastqs
            }
        }

        Array[File] fastqs_1 = select_first([fq_organize.fastqs_1, fastqs_left])
        Array[File] fastqs_2 = select_first([fq_organize.fastqs_2, fastqs_right])

        # 4. Check quality of fastq files
        call fq_qc_module.fq_qc {
            input:
                sample_id = sample_id,
                fastqs_left = fastqs_1,
                fastqs_right = fastqs_2
        }

        # 5. Align reads
        call fq_bwa_align_module.fq_bwa_align {
            input:
                fastqs_left = fastqs_1,
                fastqs_right = fastqs_2,
                sample_id = sample_id
        }

        # 6. Discarde reads with fraction of mismatches greater than or equal than mismatch_threshold
        if(add_bam_filter) {
            call bam_filter_contam_module.bam_filter_contam {
                input:
                    bam = fq_bwa_align.recalibrated_markdup_bam,
                    bai = fq_bwa_align.recalibrated_markdup_bai,
                    sample_id = sample_id,
                    mismatch_threshold = mismatch_threshold
            }
        }
    }

    # Start analysis from variant calling
    File bam_to_var_calling = select_first([bam_filter_contam.filtered_bam, fq_bwa_align.recalibrated_markdup_bam, bam])
    File bai_to_var_calling = select_first([bam_filter_contam.filtered_bai, fq_bwa_align.recalibrated_markdup_bai, bam_bai])

    # 7. Call variants
    call bam_varcalling_module.bam_varcalling {
        input:
            input_bam = bam_to_var_calling,
            input_bai = bai_to_var_calling,
            sample_id = sample_id,
            interval_list = interval_file,
            max_no_pieces_to_scatter_an_interval_file = max_no_pieces_to_scatter_an_interval_file
    }

    # 8. Filter variants
    #call vcf_var_filter_module.vcf_var_filter {
    #    input:
    #        vcf_gz = bam_varcalling.vcf_gz,
    #        vcf_gz_tbi = bam_varcalling.vcf_gz_tbi,
    #        sample_id = sample_id,
    #        interval_list = interval_file,
    #        apply_ad_filter = apply_ad_filter,
    #        ad_binom_threshold = ad_binom_threshold,
    #        analysis_type = genome_or_exome
    #}

    # 9. Annotate and filter variants
    #call vcf_anno_module.vcf_anno {
    #    input:
    #        vcf_gz = vcf_var_filter.filtered_vcf_gz,
    #        vcf_gz_tbi = vcf_var_filter.filtered_vcf_gz_tbi,
    #        vcf_anno_freq_genome_or_exome = genome_or_exome,
    #        gnomad_coverage_genome_or_exome = genome_or_exome,
    #        vcf_basename = sample_id
    #}

    # 10. Annotate variants acroding to ACMG recomendation
    # call vcf_acmg_report_module.vcf_acmg_report {
    #    input:
    #        panel_json = gene_panel_json,
    #        phenotypes_json = gene_phenotypes_json,
    #        sample_json = sample_json,
    #        patient_json = patient_json,
    #        vcf_gz = vcf_anno.annotated_and_filtered_vcf,
    #        vcf_gz_tbi = vcf_anno.annotated_and_filtered_vcf_tbi,
    #        bam = bam_to_var_calling,
    #        bai = bai_to_var_calling,
    #        realigned_bam = bam_varcalling.haplotype_caller_bam,
    #        realigned_bai = bam_varcalling.haplotype_caller_bai,
    #        other_bams = other_bams,
    #        other_bais = other_bais,
    #        sample_id = sample_id,
    #        genome_or_exome = genome_or_exome
    #}

    # 11. Call structural variants
    #if (add_snp_data_to_sv) {
    #    File? vcf_to_sv = vcf_var_filter.filtered_vcf_gz
    #    File? vcf_tbi_to_sv = vcf_var_filter.filtered_vcf_gz_tbi
    #}

    #if(genome_or_exome == "genome" && run_sv_calling) {

    #    call sv_calling_module.sv_calling {
    #        input:
    #            del_cov_max = del_cov_max,
    #            dup_cov_min = dup_cov_min,
    #            del_cov_gcbin_max = del_cov_gcbin_max,
    #            dup_cov_gcbin_min = dup_cov_gcbin_min,
    #            create_pictures = true,
    #            exclude_lcr = exclude_lcr,
    #            bam = bam_to_var_calling,
    #            bai = bai_to_var_calling,
    #            vcf_gz = vcf_to_sv,
    #            vcf_gz_tbi = vcf_tbi_to_sv,
    #            gene_panel = gene_panel_json,
    #            sample_id = sample_id
    #    }
    #}

    # 12. Estimate detection chance
    call detection_chance_module.detection_chance {
        input:
            sample_id = sample_id,
            sample_gvcf_gz = bam_varcalling.gvcf_gz,
            sample_gvcf_gz_tbi = bam_varcalling.gvcf_gz_tbi,
            panel_json = gene_panel_json,
            bam = bam_to_var_calling,
            bai = bai_to_var_calling
    }

    # 13. Check quality of bam files
    call bam_qc_module.bam_qc {
        input:
            sample_id = sample_id,
            bam = fq_bwa_align.recalibrated_markdup_bam,
            bai = fq_bwa_align.recalibrated_markdup_bai,
            genome_or_exome = genome_or_exome,
            kit = kit
    }

    # 14. Genetic sex verification
    # call sex_check_module.sex_check {
    #    input:
    #        bam = bam_to_var_calling,
    #        bai = bai_to_var_calling,
    #        vcf_gz = bam_varcalling.vcf_gz,
    #        vcf_gz_tbi = bam_varcalling.vcf_gz_tbi,
    #        sample_id = sample_id
    # }

    #  # 15. Merge pdf reports
    #  Array[File] reports_pdf = [vcf_acmg_report.pdf_report, fq_qc.quality_check_report_pdf, sex_check.report_pdf]
    #  Array[File] reports_ang_pdf = [vcf_acmg_report.ang_pdf_report, fq_qc.quality_check_report_pdf, sex_check.report_pdf]
    #
    #  call pdf_merge_task.pdf_merge {
    #    input:
    #      sample_id = sample_id,
    #      pdf = reports_pdf,
    #      ang_pdf = reports_ang_pdf
    #  }

    # 16. Merge BCO and prepare report pdf
    # Array[File] bcos_module = select_all([panel_generate.bco, panel_generate_acmg_panel.bco, resources_kit.bco, fq_organize.bco, fq_qc.bco, fq_bwa_align.bco, bam_filter_contam.bco, bam_varcalling.bco, vcf_var_filter.bco, vcf_anno.bco, vcf_acmg_report.bco, sv_calling.bco, sex_check.bco, detection_chance.bco])
    # Array[File] stdout_module = select_all([panel_generate.stdout_log, panel_generate_acmg_panel.stdout_log, resources_kit.stdout_log,  fq_organize.stdout_log, fq_qc.stdout_log, fq_bwa_align.stdout_log, bam_filter_contam.stdout_log, bam_varcalling.stdout_log, vcf_var_filter.stdout_log, vcf_anno.stdout_log, vcf_acmg_report.stdout_log, sv_calling.stdout_log, sex_check.stdout_log, detection_chance.stdout_log])
    # Array[File] stderr_module = select_all([panel_generate.stderr_log, panel_generate_acmg_panel.stderr_log, resources_kit.stderr_log, fq_organize.stderr_log, fq_qc.stderr_log, fq_bwa_align.stderr_log, bam_filter_contam.stderr_log, bam_varcalling.stderr_log, vcf_var_filter.stderr_log, vcf_anno.stderr_log, vcf_acmg_report.stderr_log, sv_calling.stderr_log, sex_check.stderr_log, detection_chance.stderr_log])

    # call bco_merge_task.bco_merge as bco_merge_pipeline {
    #    input:
    #        bco_array = bcos_module,
    #        stdout_array = stdout_module,
    #        stderr_array = stderr_module,
    #        pipeline_name = pipeline_name,
    #        pipeline_version = pipeline_version
    #}

    # call report_bco_task.report_bco as report_bco_pipeline {
    #    input:
    #        sample_id = sample_id,
    #        bco_json = bco_merge_pipeline.bco
    # }

    output {

        # 4. Check quality of fastq files
        Array[File]? fastqc_1_zips = fq_qc.fastqc_1_zips
        Array[File]? fastqc_2_zips = fq_qc.fastqc_2_zips
        File? quality_check_pdf = fq_qc.quality_check_report_pdf
        File? quality_check_html = fq_qc.quality_check_report_html

        # 5./6. Align (and filter) bam files
        File? final_bam = bam_to_var_calling
        File? final_bai = bai_to_var_calling

        # 7. Call variants
        File? gvcf_gz = bam_varcalling.gvcf_gz
        File? gvcf_gz_tbi = bam_varcalling.gvcf_gz_tbi

        File? vcf_gz = bam_varcalling.vcf_gz
        File? vcf_gz_tbi = bam_varcalling.vcf_gz_tbi

        File? final_realigned_bam = bam_varcalling.haplotype_caller_bam
        File? final_realigned_bai = bam_varcalling.haplotype_caller_bai

        # 8. Filter variants
        # File? filtered_vcf_gz = vcf_var_filter.filtered_vcf_gz
        # File? filtered_vcf_gz_tbi = vcf_var_filter.filtered_vcf_gz_tbi

        # Array[File]?  rejected_variants_vcfs = vcf_var_filter.rejected_variants_vcfs
        # Array[File]?  rejected_variants_tbis = vcf_var_filter.rejected_variants_tbis

        # File? detail_metrics_file = vcf_var_filter.detail_metrics_file
        # File? summary_metrics_file = vcf_var_filter.summary_metrics_file


        # 9. Annotate and filter variants
        # File annotated_vcf = vcf_anno.annotated_and_filtered_vcf
        # File annotated_vcf_tbi = vcf_anno.annotated_and_filtered_vcf_tbi

        # 10. Annotate variants acroding to ACMG recomendation
        # File annotated_acmg_vcf_gz = vcf_acmg_report.annotated_acmg_vcf_gz
        # File annotated_acmg_vcf_gz_tbi = vcf_acmg_report.annotated_acmg_vcf_gz_tbi
        # File? igv_screenshots_tar_gz = vcf_acmg_report.igv_screenshots_tar_gz
        # File csv_report = vcf_acmg_report.csv_report

        # File pdf_report = vcf_acmg_report.pdf_report
        # File html_report = vcf_acmg_report.html_report

        # File ang_pdf_report = vcf_acmg_report.ang_pdf_report
        # File ang_html_report = vcf_acmg_report.ang_html_report

        # 11. Call structural variants
        # File? sv_vcf = sv_calling.sv_vcf
        # File? sv_vcf_tbi = sv_calling.sv_vcf_tbi
        # File? bam_stats_json = sv_calling.bam_stats_json
        # File? sv_csv = sv_calling.csv_report

        # File? annotated_sv_vcf_gz = sv_calling.annotated_sv_vcf_gz
        # File? annotated_sv_vcf_gz_tbi = sv_calling.annotated_sv_vcf_gz_tbi
        # Array[File]? annotated_tsv_gz = sv_calling.annotated_tsv_gz
        # File? user_gene_file = sv_calling.user_gene_file

        # File? annotated_filtered_sv_vcf_gz = sv_calling.annotated_filtered_sv_vcf_gz
        # File? annotated_filtered_sv_vcf_gz_tbi = sv_calling.annotated_filtered_sv_vcf_gz_tbi

        # File? igv_html = sv_calling.igv_html
        # File? igv_pngs = sv_calling.igv_pngs

        # 12. Estimate detection chance
        File detection_chance_report_pdf = detection_chance.detection_chance_report_pdf
        File detection_chance_report_odt = detection_chance.detection_chance_report_odt
        File detection_chance_report_docx = detection_chance.detection_chance_report_docx
        File detection_chance_report_html = detection_chance.detection_chance_report_html

        # 13. Check quality of bam files
        File coverage_report_pdf = bam_qc.coverage_report_pdf
        File coverage_report_odt = bam_qc.coverage_report_odt
        File coverage_report_docx = bam_qc.coverage_report_docx
        File coverage_report_html = bam_qc.coverage_report_html

        # 14. Genetic sex verification
        # File sex_check_report_html = sex_check.report_html
        # File sex_check_report_pdf = sex_check.report_pdf

        # 15. Merge pdf reports
        #    File full_report_pdf = pdf_merge.full_report_pdf
        #    File full_report_ang_pdf = pdf_merge.full_report_ang_pdf
        #    Array[File] all_reports_pdf = pdf_merge.all_reports_pdf

        # 16. Merge BCO and prepare report pdf
        # File bco = bco_merge_pipeline.bco
        # File bco_report_html = report_bco_pipeline.bco_report_html
        # File bco_report_pdf = report_bco_pipeline.bco_report_pdf
        # File stdout_log = bco_merge_pipeline.stdout_log
        # File stderr_log = bco_merge_pipeline.stderr_log

    }
}
