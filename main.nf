#!/usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO
// Modified for chromosome-based parallelization

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

nextflow.enable.dsl = 2
//nextflow.preview.output = true

/*
========================================================================================
    PARAMETER DEFINITIONS
========================================================================================
*/


// Input/Output
params.input_folder        = null
params.input_file          = null
params.output_folder       = "pavoa_" + new Date().format('yyyyMMdd')

// Reference genome
params.ref                 = null
params.bwa_index_dir       = null
    
// Resources
params.cpu                 = 8
params.mem                 = 32
params.cpu_bqsr            = 2
params.mem_bqsr            = 10
params.cpu_dupcaller       = 64
params.mem_dupcaller       = 32
    
// Params when using folder input
params.lane                = "L001"
params.pl                  = "PL:ILLUMINA"
params.normal              = null
params.fastq_ext           = "fastq.gz"
params.suffix1             = "_1"
params.suffix2             = "_2"
    
// Known sites
params.known_sites         = []
params.mask                = null
    
// Tool options
params.bwa_option_m        = false
params.adapterremoval_opt  = ""
params.postaltjs           = null
    
//Trimming
params.trim                = false
params.adapter             = "illumina"
params.length              = 30
params.quality             = 30

// Features
params.bqsr                  = false
params.umi                   = false
params.alt                   = false

// QC
params.feature_file        = 'NO_FILE'
params.multiqc_config      = 'NO_FILE'

// Annotation
params.annovarDBlist  = null
params.annovarDBpath  = "/data/databases/annovar/hg38db/"
params.annovarBinPath = "~/bin/annovar/"
params.pass = "'PASS'"
    
// Help message
params.help                = false


/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    ========================================================================
     Mutations Pipeline v2.0
    ========================================================================
    
    Usage:
    nextflow run chromosome-parallel-alignment.nf [options]
    
    ==MANDATORY ARGUMENTS==

    Input: Choose one of the two options
    --input_file       FILE   Tab-separated file with sample info (SM, RG, pair1, pair2, normal)
    --input_folder     DIR    Folder containing FASTQ files

    Reference Genome:
    --ref              FILE   Reference genome FASTA file
    --bwa_index_dir    DIR    Directory for BWA index files (default: same as reference)

    ==OPTIONAL ARGUMENTS==

    Output:
    --output_folder   DIR     Output directory (default: bam_YYYYMMDD)
    
    Resource Options:
    --cpu             INT     CPUs for alignment (default: 8)
    --mem             INT     Memory for alignment in GB (default: 32)
    --cpu_bqsr        INT     CPUs for BQSR (default: 2)
    --mem_bqsr        INT     Memory for BQSR in GB (default: 10)
    
    Triming options:
    --trim                    Enable adapter trimming with trim_galore
    --adapter         STRING  adapter type (illumina, nextera, etc.) default(illumina)
    --length          INT      Minimum read length after trimming (default: 30)
    --quality         INT     Minimum read quality after trimming (default: 30)

    Processing Options:
    --bqsr                    Enable base quality score recalibration
    --umi            [STRING] Enable UMI-aware duplicate marking. You may pass a TAG or keep default: 'NNNNNNNN'
    --alt                     For alt aware alignment, .alt file must be present in the same directory as the reference FASTA.
    
    Quality Control:
    --feature_file    FILE    Feature file for Qualimap
    --multiqc_config  FILE    MultiQC configuration file
    
    Params when using folder input:
    pl                STRING  @RG tags, default(PL:ILLUMINA)
    lane              STRING  @RG tags, default(L001)
    normal            STRING  normal sample for somatic calling, default(null)
    fastq_ext         STRING  fastq files extentions, default(fastq.gz)
    suffix1           STRING  forward suffix, default(_1)
    suffix2           STRING  reverse suffix, default(_2)

    Known Sites and Masking:
    --known_sites      FILE   dbSNP VCF file (default: dbsnp.vcf)
    --mask             FILE   BED file for masking (default: NO_BED)

    Annotation:
    --annovarDBlist    FILE       File with two columns : protocols and operations (see anovar documentation).'
    --annovarDBpath    PATH       Path to annovarDB.'
    --annovarBinPath   PATH       Path to table_annovar.pl.'
    --pass             STRING     filter flags, comma separated list'

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    PARAMETER VALIDATION
========================================================================================
*/

// Validate required parameters
if (!params.input_folder && !params.input_file) {
    error "Either --input_folder or --input_file must be specified"
}

if (!params.ref) {
    error "--ref must be specified"
}

outputDir = "${params.output_folder}"

/*
========================================================================================
    LOGGING 
========================================================================================
*/

// Log parameters
log.info """
========================================================================
 Alignment Pipeline v2.0
========================================================================
Input folder       : ${params.input_folder}
Input file         : ${params.input_file}
Reference file     : ${params.ref}
bwa_index_dir      : ${params.bwa_index_dir ?: 'Same as reference'}
Output folder      : ${params.output_folder}
CPUs               : ${params.cpu}
CPUs for BQSR      : ${params.cpu_bqsr}
Memory             : ${params.mem} GB
Memory for BQSR    : ${params.mem_bqsr} GB
Recalibration      : ${params.bqsr}
Trimming           : ${params.trim}
Adapter            : ${params.adapter}
Length             : ${params.length}
Quality            : ${params.quality}
UMI                : ${params.umi}
Known sites        : ${params.known_sites}
Feature file       : ${params.feature_file}
MultiQC config     : ${params.multiqc_config}
Plateform          : ${params.pl}
fastq_ext          : ${params.fastq_ext}
suffix1            : ${params.suffix1}
suffix2            : ${params.suffix2}
annovarDBlist      : ${params.annovarDBlist}
annovarDBpath      : ${params.annovarDBpath}
annovarBinPath     : ${params.annovarBinPath}
pass               : ${params.pass}

========================================================================
    Workflow started at: ${new Date().format('yyyy-MM-dd HH:mm:ss')}
========================================================================
"""
    


/*
========================================================================================
    INCLUDE MODULES
========================================================================================
*/

include { PREPARE_REFERENCES           } from './modules/references'
include { get_chromosomes              } from './modules/references'
include { PREPARE_REFERENCES_DUPCALLER } from './modules/references'
include { PREPARE_INPUT_FROM_TSV       } from './modules/input'
include { PREPARE_INPUT_FROM_FOLDER    } from './modules/input'
include { PREPARE_CALLING_INPUT        } from './modules/input'
include { DUPCALLER_TRIM               } from './modules/dupcaller'
include { DUPCALLER_CALL               } from './modules/dupcaller'
include { TRIM_GALORE                  } from './modules/trimming'
include { FASTQ_ALIGNMENT              } from './modules/alignment'
include { MARK_DUPLICATES              } from './modules/gatk'
include { BQSR                         } from './modules/gatk'
include { MERGE_BAM                    } from './modules/merging'
include { QUALIMAP                     } from './modules/qc'
include { FLAGSTAT                     } from './modules/qc'
include { MULTIQC                      } from './modules/qc'
include { ANNOTATION                   } from './modules/annovar'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    main:

    // Define input channels
    samples = Channel.empty() // Initialize empty channel for reports
    reports = Channel.empty() // Initialize empty channel for reports

    // Prepare input samples
    if (params.input_file) {
        PREPARE_INPUT_FROM_TSV(params.input_file)
        samples = PREPARE_INPUT_FROM_TSV.out.reads
    } else {
        PREPARE_INPUT_FROM_FOLDER(params.input_folder, params.fastq_ext, params.suffix1, params.suffix2)
        samples = PREPARE_INPUT_FROM_FOLDER.out.reads
    }

    // Prepare chromosome-specific reference files
    PREPARE_REFERENCES(params.ref, params.bwa_index_dir, params.output_folder)
    ref = PREPARE_REFERENCES.out.ref
    indexes = PREPARE_REFERENCES.out.indexes

    // Prepare vcf files for known sites (dbsnp, dbindel, etc.)
    known_sites = params.known_sites ? Channel
        .fromList(params.known_sites)
        .map { vcf -> [file(vcf), file("${vcf}.tbi")] }
        .collect()
    : Channel.empty()

    // Prepare mask and feature file
    mask = params.mask ? (file(params.mask)) : (file("NO_BED"))
    feature_file = (params.feature_file) ? file(params.feature_file) : file("NO_FEATURE")

    // DUPCALLER TRIMMING if UMI is enabled
    if (params.umi) {
        DUPCALLER_TRIM(samples)
        samples = DUPCALLER_TRIM.out.trimmed_fastq
    }

    // Trim_galore if enabled
    if (params.trim) {
        TRIM_GALORE(samples)
        samples = TRIM_GALORE.out.trimmed_fastq
        reports = reports.concat(TRIM_GALORE.out.reports)
    }

    // Align FASTQ files to reference genome
    FASTQ_ALIGNMENT(samples, ref, indexes)

    // Merge chromosome BAMs per sample
    alignments_ch = FASTQ_ALIGNMENT.out.bam_files
        .map { sample_id, rg, bams, bais -> [sample_id, bams, bais] }
        .groupTuple(by: 0) // Group by sample ID

    MERGE_BAM(alignments_ch)
    alignments_merged_ch = MERGE_BAM.out.merged_bam

    // Mark duplicates per chromosome (if UMI enabled)
    if (params.umi) {
        MARK_DUPLICATES(alignments_merged_ch)
        alignments_merged_ch = MARK_DUPLICATES.out.bam_files
    }
    
    // Base quality score recalibration per chromosome
    if (params.bqsr) {
        BQSR(alignments_merged_ch, ref, indexes, known_sites)
        alignments_merged_ch = BQSR.out.bamfiles
        reports = reports.concat(BQSR.out.recal_table_files)
    }
    
    // Qualimap for BAM files
    QUALIMAP(alignments_merged_ch, feature_file)
    reports = reports.concat(QUALIMAP.out.report)

    // Flagstat for BAM files
    FLAGSTAT(alignments_merged_ch)
    reports = reports.concat(FLAGSTAT.out.report)

    // MultiQC reports
    MULTIQC(reports.collect(), file(params.multiqc_config))

    // Stop here if input_file_path is null
    //if (!params.input_file) { exit() }

    // Organise Samples for Calling
    PREPARE_CALLING_INPUT(params.input_file, alignments_merged_ch)
    pairs = PREPARE_CALLING_INPUT.out.pairs

    // pairs.view()
    // Call variants using DUPCALLER
    chromosome = get_chromosomes(ref,indexes).splitText().collect{ it.trim() }
    DUPCALLER_CALL(pairs, ref, indexes, chromosome, known_sites, mask)
    vcfs = DUPCALLER_CALL.out.vcfs

    // Annotate VCF files
    if(params.annovarDBlist){
        ANNOTATION(vcfs)
    }

    //publish:
    //final_alignment = alignments_merged_ch
    publish_final_alignment(alignments_merged_ch)


}

process publish_final_alignment {

    tag "${sample_id}"
    label 'publish'

    publishDir "${params.output_folder}/BAM/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_files), path(bai_files)

    output:
    tuple val(sample_id), path(bam_files), path(bai_files)

    script:
    """
    """
}   

/*output {

    final_alignment {
        path 'BAM'
        mode 'copy'
    }
}*/

/*
========================================================================================
    RUN only dupcaller
========================================================================================
*/

workflow dupcaller{

    main:

    // Define input channels
    pairs = Channel
        .fromPath(params.input_file)
        .splitCsv(header: true, sep: '\t', strip: true)
        .filter { row -> row.normal } // Only keep samples with matched normals
        .map { row -> tuple(row.SM, file(row.tumor), file(row.tumor + ".bai"), file(row.normal), file(row.normal + ".bai")) }

    // Prepare chromosome-specific reference files
    PREPARE_REFERENCES_DUPCALLER(params.ref, params.output_folder)
    ref = PREPARE_REFERENCES_DUPCALLER.out.ref
    indexes = PREPARE_REFERENCES_DUPCALLER.out.indexes

    // Prepare vcf files for known sites (dbsnp, dbindel, etc.)
    known_sites = params.known_sites ? Channel
        .fromList(params.known_sites)
        .map { vcf -> [file(vcf), file("${vcf}.tbi")] }
        .collect()
    : Channel.empty()

    // Prepare mask and feature file
    mask = params.mask ? (file(params.mask)) : (file("NO_BED"))

    // pairs.view()
    // Call variants using DUPCALLER
    DUPCALLER_CALL(pairs, ref, indexes, known_sites, mask)
    vcfs = DUPCALLER_CALL.out.vcfs

    // Annotate VCF files
    if(params.annovarDBlist){
        ANNOTATION(vcfs)
    }

}




/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    ========================================================================
     Pipeline completed successfully!
    ========================================================================
    Results directory: ${params.output_folder}
    Duration: ${workflow.duration}
    Success: ${workflow.success}
    Exit status: ${workflow.exitStatus}
    ========================================================================
    """.stripIndent()
}
