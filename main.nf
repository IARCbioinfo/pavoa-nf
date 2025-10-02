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
params.input_file          = null
params.input_folder        = null  // not used with current version
params.output_folder       = "pavoa_output"

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

// Params when using folder input (not used with current version)
params.lane                = "L001"
params.normal              = null
params.fastq_ext           = "fastq.gz"
params.suffix1             = "_1"
params.suffix2             = "_2"

//Trimming options
//params.adapterremoval_opt  = ""
params.umi                 = false
params.trim                = true
params.adapter             = "illumina"
params.length              = 30
params.quality             = 30


//Mapping options
params.bwa_option_m        = true
params.pl                  = "PL:ILLUMINA"
params.postaltjs           = null
params.bqsr                = params.umi ? false : true
params.known_sites         = []
params.snp_contam          = null
params.recall              = false  // Skip alignment step and use existing BAM files in project folder.

// QC
params.feature_file        = 'NO_FILE'
params.multiqc_config      = 'NO_FILE'

//dupcaller 
params.mask                = null
// params.chromosome = null

// strelka2
params.strelka2 = true
params.strelka_bin = "/opt/conda/envs/strelka2-nf/share/strelka-2.9.10-0/bin/"
params.strelka_config = null
params.exome = false

// mutect2
params.mutect2 = true
params.mutect_args = ""
params.PON = null
params.nsplit = 1000

// Annotation
params.annovarDBlist  = null
params.annovarDBpath  = "/data/databases/annovar/hg38db/"
params.annovarBinPath = "~/bin/annovar/"
params.pass = "'PASS'"

// Filtering parameters
params.cov_n_thresh        = params.umi ? 1 : 10      // Normal coverage threshold
params.cov_t_thresh        = params.umi ? 1 : 10      // Tumor coverage threshold  
params.min_vaf_t_thresh    = params.umi ? 0 : 0.1     // Min Tumor VAF threshold
params.max_vaf_t_thresh    = 1                        // Max Tumor VAF threshold
params.cov_alt_t_thresh    = params.umi ? 1 : 3       // Tumor alternative allele coverage threshold

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

    Input:
    --input_file      FILE    Tab-separated file with sample info (SM, RG, pair1, pair2, normal)

    Reference Genome:
    --ref             FILE    Reference genome FASTA file
    --bwa_index_dir   DIR     Directory for BWA index files (default: same as reference)

    ==OPTIONAL ARGUMENTS==

    Output:
    --output_folder   DIR     Output directory (default: pavoa_output)
    
    Resource Options:
    --cpu             INT     CPUs for alignment (default: 8)
    --mem             INT     Memory for alignment in GB (default: 32)
    --cpu_bqsr        INT     CPUs for BQSR (default: 2)
    --mem_bqsr        INT     Memory for BQSR in GB (default: 10)
    --cpu_dupcaller   INT     CPUs for DupCaller (default: 64)
    --mem_dupcaller   INT     Memory for DupCaller in GB (default: 32)

    Triming options:
    --umi             STRING  Enable UMI-aware duplicate marking. You may pass a TAG or keep default: 'NNNNNNNN'
    --trim                    Enable adapter trimming with trim_galore (default: true)
    --adapter         STRING  adapter type (illumina, nextera, etc.) default(illumina)
    --length          INT     Minimum read length after trimming (default: 30)
    --quality         INT     Minimum read quality after trimming (default: 30)

    Mapping options:
    --bwa_option_m            Use -M option in BWA and Samblaster (for Picard compatibility)
    --pl              STRING  Platform for RG group (default: "PL:ILLUMINA")
    --postaltjs       FILE    Path to bwa-postalt.js script (default: /opt/bwa-postalt.js)
    --bqsr                    Enable base quality score recalibration (default: disabled if UMI is enabled)
    --known_sites     FILE    dbSNP VCF file, also use for calling (default: dbsnp.vcf)
    --snp_contam      FILE    SNP VCF file for ContEst (default: null)
    --recall                  Skip alignment step and use existing BAM files in project folder.
    
    Quality Control:
    --feature_file    FILE    Feature file for Qualimap
    --multiqc_config  FILE    MultiQC configuration file

    Dupcaller:
    --mask             FILE   BED file for masking (default: NO_BED)

    Strelka2:
    --strelka2                Enable Strelka2 caller (default: true)
    --strelka_bin     DIR     Path to Strelka2 bin directory (default: /opt/conda/envs/strelka2-nf/share/strelka-2.9.10-0/bin/)
    --strelka_config  FILE    Path to Strelka2 config file (default: null)
    --exome                   Use exome mode for Strelka2 (default: false)

    Mutect2:
    --mutect2                 Enable Mutect2 caller (default: true)
    --mutect_args    STRING   Additional arguments for Mutect2 (default: "")
    --nsplit         INT      Number of chromosomes to process in parallel (default: 1000)

    Annotation:
    --annovarDBlist    FILE       File with two columns : protocols and operations (see anovar documentation).'
    --annovarDBpath    PATH       Path to annovarDB.'
    --annovarBinPath   PATH       Path to table_annovar.pl.'
    --pass             STRING     filter flags, comma separated list, default ("PASS")'

    Filtering :
    --cov_n_thresh     INT    Coverage threshold in Normal sample (default: 10)
    --cov_t_thresh     INT    Coverage threshold in Tumor sample(default: 10)
    --min_vaf_t_thresh FLOAT  Min VAF threshold in Tumor sample (default: 0.1)
    --max_vaf_t_thresh FLOAT  Max VAF threshold in Tumor sample (default: 1)
    --cov_alt_t_thresh INT    Alternative allele coverage threshold in sample (default: 3)

    """.stripIndent()
}




/*
========================================================================================
    LOG PARAMETERS VALUES 
========================================================================================
*/
def logParameters() {

    // Log parameters
    log.info """
    ========================================================================
    Alignment Pipeline v2.0
    ========================================================================
    Input file           : ${params.input_file}
    Output folder        : ${params.output_folder}
    Reference file       : ${params.ref}

    CPUs                 : ${params.cpu}
    CPUs for BQSR        : ${params.cpu_bqsr}
    CPUs for DupCaller   : ${params.cpu_dupcaller}
    Memory               : ${params.mem} GB
    Memory for BQSR      : ${params.mem_bqsr} GB
    Memory for DupCaller : ${params.mem_dupcaller} GB

    UMI                  : ${params.umi}
    Trimming             : ${params.trim}
    Adapter              : ${params.adapter}
    Length               : ${params.length}
    Quality              : ${params.quality}

    BWA -M option        : ${params.bwa_option_m}
    Plateform            : ${params.pl}
    Postaltjs            : ${params.postaltjs}
    BQSR                 : ${params.bqsr}
    Known sites          : ${params.known_sites}
    SNP contam           : ${params.snp_contam}
    reCalling only       : ${params.recall}

    Feature file         : ${params.feature_file}
    MultiQC config       : ${params.multiqc_config}

    Mask                 : ${params.mask}
    Strelka2             : ${params.strelka2}
    Strelka bin          : ${params.strelka_bin}
    Strelka config       : ${params.strelka_config}
    Exome mode           : ${params.exome}
    Mutect2              : ${params.mutect2}
    Mutect2 args         : ${params.mutect_args}
    PON                  : ${params.PON}
    Nsplit               : ${params.nsplit}

    annovarDBlist        : ${params.annovarDBlist}
    annovarDBpath        : ${params.annovarDBpath}
    annovarBinPath       : ${params.annovarBinPath}
    pass                 : ${params.pass}

    Normal cov thresh    : ${params.cov_n_thresh}
    Tumor cov thresh     : ${params.cov_t_thresh}
    Min VAF thresh       : ${params.min_vaf_t_thresh}
    Max VAF thresh       : ${params.max_vaf_t_thresh}
    alt cov thresh       : ${params.cov_alt_t_thresh}
    ========================================================================
    """
}


/*
========================================================================================
    INCLUDE MODULES
========================================================================================
*/

include { PREPARE_REFERENCES           } from './modules/references'
include { PREPARE_REFERENCES_DUPCALLER } from './modules/references'
include { PREPARE_INPUT_FROM_TSV       } from './modules/input'
include { PREPARE_INPUT_FROM_FOLDER    } from './modules/input'
include { PREPARE_CALLING_INPUT        } from './modules/input'
include { DUPCALLER_TRIM               } from './modules/dupcaller'
include { DUPCALLER_CALL               } from './modules/dupcaller'
include { DUPCALLER_ESTIMATE           } from './modules/dupcaller'
include { TRIM_GALORE                  } from './modules/trimming'
include { FASTQ_ALIGNMENT              } from './modules/alignment'
include { MARK_DUPLICATES              } from './modules/gatk'
include { BQSR                         } from './modules/gatk'
include { MERGE_BAM                    } from './modules/merging'
include { QUALIMAP                     } from './modules/qc'
include { FLAGSTAT                     } from './modules/qc'
include { MULTIQC                      } from './modules/qc'
include { STRELKA2_CALL                } from './modules/strelka'
include { MUTECT2_CALL                 } from './modules/mutect2'
include { ANNOTATION                   } from './modules/annovar'
include { GAMA_FILTER                  } from './modules/annovar'



/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    main:

    if (params.help) {
        helpMessage()
        exit 0
    }

    //PARAMETER VALIDATION
    logParameters()

    //check if params.input_file exists
    if (!params.input_file) {
        error "Either --input_folder or --input_file must be specified"
        def input_file = file(params.input_file)
        if (!input_file.exists()) {
            error "Input file ${params.input_file} does not exist"  
        }
    }

    //check if params.ref is specified and exists
    if (!params.ref) {
        error "--ref must be specified"
        def ref_file = file(params.ref)
        if (!ref_file.exists()) {
            error "Reference file ${params.ref} does not exist"
        }
    }

    //check if params.recall true, then params.output_folder must exist and contain BAM files
    if (params.recall) {
        def bam_files = file("${params.output_folder}/BAM/*.bam")
        if (bam_files.size()==0) {
            error "When --recall is enabled, the output folder must contain existing BAM files in the BAM subfolder."
        }
    }

    //check if known_sites files exist
    if (params.known_sites) {
        params.known_sites.each { vcf ->
            def vcf_file = file(vcf)
            if (!vcf_file.exists()) {
                error "Known sites file ${vcf} does not exist"
            }
            def tbi_file = file("${vcf}.tbi")
            if (!tbi_file.exists()) {
                error "Index file for known sites ${vcf}.tbi does not exist"
            }
        }
    }

    //check if feature_file exists
    if (params.feature_file != 'NO_FILE') {
        def feature_file = file(params.feature_file)
        if (!feature_file.exists()) {
            error "Feature file ${params.feature_file} does not exist"
        }
    }

    // Check if all parameters for annovar are specified if annovarDBlist is specified
    if (params.annovarDBlist) {
        if (!params.annovarDBpath) {
            error "--annovarDBpath must be specified when --annovarDBlist is specified"
        }
        if (!params.annovarBinPath) {
            error "--annovarBinPath must be specified when --annovarDBlist is specified"
        }
        def annovarDBlist_file = file(params.annovarDBlist)
        if (!annovarDBlist_file.exists()) {
            error "Annovar DB list file ${params.annovarDBlist} does not exist"
        }
    }   

    // INPUT CHANNELS
    samples = Channel.empty() // Initialize empty channel for reports
    
    // Prepare input samples
    PREPARE_INPUT_FROM_TSV(params.input_file)
    samples = PREPARE_INPUT_FROM_TSV.out.reads

    // Prepare reference files
    PREPARE_REFERENCES(params.ref, params.bwa_index_dir, params.output_folder)
    ref = PREPARE_REFERENCES.out.ref
    indexes = PREPARE_REFERENCES.out.indexes

    // Prepare vcf files for known sites (dbsnp, dbindel, etc.)
    known_sites = params.known_sites ? Channel
        .fromList(params.known_sites)
        .map { vcf -> [file(vcf), file("${vcf}.tbi")] }
        .collect()
    : Channel.empty()

    // Prepare feature file
    feature_file = (params.feature_file) ? file(params.feature_file) : file("NO_FEATURE")

    if(!params.recall){
        // MAPPING
        MAPPING(samples, ref, indexes, known_sites, feature_file)
        alignments = MAPPING.out.final_alignment.collect()
    }else{
        // Instead of running MAPPING, collect existing BAM files
        alignments = channel.fromPath( params.output_folder + "/BAM/*.bam" )
            .map { bam -> 
            def bam_name = bam.baseName
                .replaceAll("_mkdup","")
                .replaceAll("_bqsr","")
            tuple(bam_name, bam, file("${bam}.bai")) }
    }

    // Organise Samples for Calling
    PREPARE_CALLING_INPUT(params.input_file, alignments)
    pairs = PREPARE_CALLING_INPUT.out.pairs

    // CALLING
    CALLING(pairs, ref, indexes, known_sites, feature_file)

}

/*
========================================================================================
    MAPPING WORKFLOW
========================================================================================
*/
workflow MAPPING {

    take:
    samples
    ref
    indexes
    known_sites
    feature_file

    main:

    reports = Channel.empty() // Initialize empty channel for reports

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

    // Group BAMs per sample
    alignments_ch = FASTQ_ALIGNMENT.out.bam_files
        .map { sample_id, _rg, bams, bais -> [sample_id, bams, bais] }
        .groupTuple(by: 0) // Group by sample ID

    MERGE_BAM(alignments_ch)
    alignments_merged_ch = MERGE_BAM.out.merged_bam

    MARK_DUPLICATES(alignments_merged_ch)
    alignments_merged_ch = MARK_DUPLICATES.out.bam_files
    
    // Base quality score recalibration
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
    reports = reports.collect()
        .map { files -> files.reverse().unique { it.getName() } } // Keep only one file per type
    MULTIQC(reports, file(params.multiqc_config))

    //publish:
    publish_final_alignment(alignments_merged_ch)

    emit:
    final_alignment = alignments_merged_ch

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
    CALLING WORKFLOW
========================================================================================
*/
workflow CALLING {

    take:
    pairs
    ref
    indexes
    known_sites
    bed

    main:

    // Collect all VCF files from different callers into one channel
    all_vcfs = Channel.empty()

    // Run DUPCALLER if UMI mode is enabled
    if (params.umi) {
        mask = params.mask ? (file(params.mask)) : (file("NO_BED"))
        DUPCALLER_CALL(pairs, ref, indexes, known_sites, mask)
        all_vcfs = all_vcfs.mix(DUPCALLER_CALL.out.vcfs)
    }

    // Run Strelka2 if enabled
    if(params.strelka2) {
        STRELKA2_CALL(pairs, ref, indexes, bed)
        all_vcfs = all_vcfs.mix(STRELKA2_CALL.out.vcfs)
    }

    // Run Mutect2 if enabled
    if(params.mutect2) {
        MUTECT2_CALL(pairs, ref, indexes, bed, known_sites)
        all_vcfs = all_vcfs.mix(MUTECT2_CALL.out.vcfs)
    }

    // Run annotation once with all VCFs
    if(params.annovarDBlist) {
        filtered_vcf = ANNOTATION(all_vcfs)

        if(params.umi) {
            // Filter for dupcaller VCFs using caller tag
            dupcaller_vcfs = filtered_vcf.filter { it[2] == 'dupcaller' }
                .map { sample_id, files, _caller_tag ->tuple(sample_id, files)}
            dupcaller_vcfs.view()

            // Prepare input for DUPCALLER_ESTIMATE
            sit = dupcaller_vcfs
                .concat(DUPCALLER_CALL.out.trinuc)
                .groupTuple(by: 0)
                .map { sample_id, files ->
                    def snv_vcf = files.find { it.name.contains('_snv') }
                    def indel_vcf = files.find { it.name.contains('_indel') }
                    def trinuc_file = files.find { it.name.contains('trinuc') }
                    tuple(sample_id, snv_vcf, indel_vcf, trinuc_file)
                }//.filter { it[2] != null }
                
            sit.view()
                
            // Reestimate duplication rates
            DUPCALLER_ESTIMATE(sit, ref, indexes)
        }

    }

}


workflow.onComplete {
    // THE END
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