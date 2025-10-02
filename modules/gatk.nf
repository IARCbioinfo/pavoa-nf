/*
========================================================================================
    POSTPROCESSING MODULE (BQSR, DUPLICATE MARKING, ETC.)
========================================================================================
*/

process MARK_DUPLICATES_UMI {
    tag "${file_tag}"
    label 'gatk'

    publishDir "${params.output_folder}/QC/duplicates/", mode: 'copy', pattern: "*_mkdup.metrics"

    cpus 2
    memory "16.GB"
    
    input:
    tuple val(file_tag), path(bam), path(bai)

    output:
    tuple val(file_tag), path("*_mkdup.bam"), path("*_mkdup.bam.bai"), emit: bam_files
    path("*_mkdup.metrics"), emit: metrics

    script:
    def file_tag_new = bam.baseName + '_mkdup'
    """
    set -euo pipefail
    
    mkdir -p tmpdir
    
    # Mark duplicates with UMI support
    gatk --java-options '-Xmx14G -Xms14G -Xss512k -XX:ReservedCodeCacheSize=512M -Djava.io.tmpdir=tmpdir' \\
        MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${file_tag_new}.tmp.bam \\
        --METRICS_FILE ${file_tag_new}_mkdup.metrics \\
        --DUPLEX_UMI true \\
        --ASSUME_SORT_ORDER coordinate \\
        --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    
    # Fix header if needed
    samtools view -H ${file_tag_new}.tmp.bam > header.sam
    sed -i 's/SO:unknown/SO:coordinate/' header.sam
    samtools reheader header.sam ${file_tag_new}.tmp.bam > ${file_tag_new}.bam
    samtools index ${file_tag_new}.bam ${file_tag_new}.bam.bai

    # Cleanup
    rm -f ${file_tag_new}.tmp.bam ${file_tag_new}.tmp.bai header.sam
    rm -rf tmpdir
    """
    
    stub:
    def file_tag_new = bam.baseName + '_mkdup'
    """
    touch ${file_tag_new}.bam ${file_tag_new}.bam.bai ${file_tag_new}.metrics
    """
}

process MARK_DUPLICATES_STANDARD {
    tag "${file_tag}"
    label 'gatk'
    
    cpus 2
    memory "16.GB"
    
    input:
    tuple val(file_tag), path(bam), path(bai)
    
    output:
    tuple val(file_tag), path("*_mkdup.bam"), path("*_mkdup.bam.bai"), emit: bam_files
    path("*_mkdup.metrics"), emit: metrics

    script:
    def file_tag_new = bam.baseName + '_mkdup'
    """
    set -euo pipefail
    
    mkdir -p tmpdir
    
    # Standard duplicate marking
    gatk --java-options '-Xmx6G -Xms6G -Djava.io.tmpdir=tmpdir' \\
        MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${file_tag_new}.bam \\
        --METRICS_FILE ${file_tag_new}.metrics \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT

    samtools index ${file_tag_new}.bam ${file_tag_new}.bam.bai

    # Cleanup
    rm -rf tmpdir
    """
    
    stub:
    def file_tag_new = bam.baseName + '_mkdup'
    """
    touch ${file_tag_new}.bam ${file_tag_new}.bam.bai ${file_tag_new}.metrics
    """
}


/***************************************************************************************/
/************************  Process : base_quality_score_recalibration *************/
/***************************************************************************************/

process BQSR {
   tag "${file_tag}"
   label 'gatk'

   cpus params.cpu_bqsr
   memory params.mem_bqsr+'GB'

   publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*bam*"
   publishDir "$params.output_folder/QC/BAM/BQSR/", mode: 'copy',
   saveAs: {filename ->
       if (filename.indexOf("table") > 0) "$filename"
       else if (filename.indexOf("plots") > 0) "$filename"
       else null
	}

  input:
    tuple val(file_tag), path(bam), path(bai)
    path(ref)
    path(indexes)
    path(known_sites)
    
  output:
    tuple val(file_tag), path("*_bqsr.bam"), path("*_bqsr.bam.bai"), emit: bamfiles
    path("*_recal.table"), emit: recal_table_files
    path("*plots.pdf"), emit: plots

  script:
    def file_name = bam.baseName
    def file_tag_new = file_name+'_bqsr'
    def known_sites_args = known_sites.findAll { it.name.endsWith('.vcf.gz') }.collect { "--known-sites ${it}" }.join(' ')
    """
    gatk BaseRecalibrator --java-options "-Xmx${task.memory.toGiga()}G" -R $ref -I $bam ${known_sites_args} -O ${file_name}_recal.table
    gatk ApplyBQSR --java-options "-Xmx${task.memory.toGiga()}G" -R $ref -I $bam --bqsr-recal-file ${file_name}_recal.table -O ${file_tag_new}.bam
    gatk BaseRecalibrator --java-options "-Xmx${task.memory.toGiga()}G" -R $ref -I ${file_tag_new}.bam ${known_sites_args} -O ${file_tag_new}_recal.table
    #TODO: check R environment in the container for this to work
    #gatk AnalyzeCovariates --java-options "-Xmx${task.memory.toGiga()}G" -before ${file_name}_recal.table -after ${file_tag_new}_recal.table -plots ${file_tag_new}_recalibration_plots.pdf
    touch ${file_tag_new}_recalibration_plots.pdf	
    mv ${file_tag_new}.bai ${file_tag_new}.bam.bai
    """

  stub:
    def file_name = bam.baseName
    def file_tag_new = file_name+'_bqsr'
    """
    touch ${file_tag_new}.bam ${file_tag_new}.bam.bai
    touch ${file_tag_new}_recalibration_plots.pdf
    touch ${file_tag_new}_recal.table
    """
}



workflow MARK_DUPLICATES {
    take:
    bam_files
    
    main:
    if (params.umi) {
        MARK_DUPLICATES_UMI(bam_files)
        marked_bams = MARK_DUPLICATES_UMI.out.bam_files
        metrics = MARK_DUPLICATES_UMI.out.metrics
    } else {
        MARK_DUPLICATES_STANDARD(bam_files)
        marked_bams = MARK_DUPLICATES_STANDARD.out.bam_files
        metrics = MARK_DUPLICATES_STANDARD.out.metrics
    }
    
    emit:
    bam_files = marked_bams
    metrics = metrics
}


