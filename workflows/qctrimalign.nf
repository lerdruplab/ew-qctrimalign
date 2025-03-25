/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAT_FASTQ as CAT_FASTQ   } from '../modules/nf-core/cat/fastq/main'
include { FASTQC as FASTQC_RAW     } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../modules/nf-core/fastqc/main'
include { MULTIQC                  } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE               } from '../modules/nf-core/trimgalore/main'
include { BOWTIE_ALIGN             } from '../modules/nf-core/bowtie/align/main'
include { BOWTIE2_ALIGN            } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORT            } from '../modules/nf-core/samtools/sort/main'
include { BEDTOOLS_BAMTOBED        } from '../modules/nf-core/bedtools/bamtobed/main'
include { paramsSummaryMap         } from 'plugin/nf-validation'
include { paramsSummaryMultiqc     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText   } from '../subworkflows/local/utils_nfcore_qctrimalign_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QCTRIMALIGN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_index       // channel: prebuilt index read in from --index
    aligner        // string: defines the aligner read in from --aligner

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
            .set { ch_fastq }

    //
    // MODULE: Cat fastq files if needed
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_fastq_filtered }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC_RAW (
        ch_fastq_filtered
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    //
    // MODULE: Run cutadapt
    //
    TRIMGALORE (
        ch_fastq_filtered
    )
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]})
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC_TRIMMED (
        TRIMGALORE.out.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())

    if (params.aligner == 'bowtie1') {

      //
      // MODULE: Run bowtie/align
      //
      BOWTIE_ALIGN (
          TRIMGALORE.out.reads,
          ch_index
      )
      ch_multiqc_files = ch_multiqc_files.mix(BOWTIE_ALIGN.out.log.collect{it[1]})
      ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)

      //
      // MODULE: Run samtools/sort
      //
      SAMTOOLS_SORT (
          BOWTIE_ALIGN.out.bam
      )
      ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    } else if (params.aligner == 'bowtie2') {

      //
      // MODULE: Run bowtie2/align
      //
      BOWTIE2_ALIGN (
          TRIMGALORE.out.reads,
          ch_index,
          false,
          false
      )
      ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]})
      ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

      //
      // MODULE: Run samtools/sort
      //
      SAMTOOLS_SORT (
          BOWTIE2_ALIGN.out.bam
      )
      ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    } else {

          error "Invalid tool choice: ${params.aligner}"

    }

    //
    // MODULE: Run samtools/sort
    //
    BEDTOOLS_BAMTOBED (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
