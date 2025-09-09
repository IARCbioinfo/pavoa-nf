/*
========================================================================================
    POSTPROCESSING MODULE (BQSR, DUPLICATE MARKING, ETC.)
========================================================================================
*/

process MARK_DUPLICATES_UMI {
    tag "${file_tag}"
    label 'gatk'

    publishDir "${params.output_folder}/QC/duplicates/", mode: 'copy', pattern: "${file_tag}_dedup.metrics"

    cpus 2
    memory "16.GB"
    
    input:
    tuple val(file_tag), path(bam), path(bai)

    output:
    tuple val(file_tag), path("${file_tag}_dedup.bam"), path("${file_tag}_dedup.bam.bai"), emit: bam_files
    path("${file_tag}_dedup.metrics"), emit: metrics

    script:
    """
    set -euo pipefail
    
    mkdir -p tmpdir
    
    # Mark duplicates with UMI support
    gatk --java-options '-Xmx14G -Xms14G -Xss512k -XX:ReservedCodeCacheSize=512M -Djava.io.tmpdir=tmpdir' \\
        MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${file_tag}_dedup.tmp.bam \\
        --METRICS_FILE ${file_tag}_dedup.metrics \\
        --DUPLEX_UMI true \\
        --ASSUME_SORT_ORDER coordinate \\
        --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    
    # Fix header if needed
    samtools view -H ${file_tag}_dedup.tmp.bam > header.sam
    sed -i 's/SO:unknown/SO:coordinate/' header.sam
    samtools reheader header.sam ${file_tag}_dedup.tmp.bam > ${file_tag}_dedup.bam
    samtools index ${file_tag}_dedup.bam -o ${file_tag}_dedup.bam.bai

    # Cleanup
    rm -f ${file_tag}_dedup.tmp.bam ${file_tag}_dedup.tmp.bai header.sam
    rm -rf tmpdir
    """
    
    stub:
    """
    touch ${file_tag}_dedup.bam ${file_tag}_dedup.bam.bai ${file_tag}_dedup.metrics
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
    tuple val(file_tag), path("${file_tag}_dedup.bam"), path("${file_tag}_dedup.bam.bai"), emit: bam_files
    path("${file_tag}_dedup.metrics"), emit: metrics

    script:
    """
    set -euo pipefail
    
    mkdir -p tmpdir
    
    # Standard duplicate marking
    gatk --java-options '-Xmx6G -Xms6G -Djava.io.tmpdir=tmpdir' \\
        MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${file_tag}_dedup.bam \\
        --METRICS_FILE ${file_tag}_dedup.metrics \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    
    mv "${file_tag}_dedup.bai" "${file_tag}_dedup.bam.bai"
    
    # Cleanup
    rm -rf tmpdir
    """
    
    stub:
    """
    touch ${file_tag}_dedup.bam ${file_tag}_dedup.bam.bai ${file_tag}_dedup.metrics
    """
}

process MARK_DUPLICATES_STANDARD_CHR {
    tag "${file_tag}"
    label 'gatk'
    
    cpus 2
    memory "8.GB"
    
    input:
    tuple val(file_tag), val(chromosome), path(bam), path(bai)
    
    output:
    tuple val(file_tag), val(chromosome), path("${output_prefix}_dedup.bam"), path("${output_prefix}_dedup.bam.bai"), emit: bam_files
    path("${output_prefix}_dedup.metrics"), emit: metrics
    
    script:
    output_prefix = "${file_tag}.${chromosome}"
    
    """
    set -euo pipefail
    
    mkdir -p tmpdir
    
    # Standard duplicate marking by chromosome
    gatk --java-options '-Xmx6G -Xms6G -Djava.io.tmpdir=tmpdir' \\
        MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${output_prefix}_dedup.bam \\
        --METRICS_FILE ${output_prefix}_dedup.metrics \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    
    # Cleanup
    rm -rf tmpdir
    """
    
    stub:
    output_prefix = "${file_tag}.${chromosome}"
    """
    touch ${output_prefix}_dedup.bam ${output_prefix}_dedup.bam.bai ${output_prefix}_dedup.metrics
    """
}

/***************************************************************************************/
/************************  Process : mebase_quality_score_recalibrationrge *************/
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
    tuple val(file_tag), path("*_BQSRecalibrated.bam"), path("*_BQSRecalibrated.bam.bai"), emit: bamfiles
    path("*_recal.table"), emit: recal_table_files
    path("*plots.pdf"), emit: plots

  shell:
    def file_name = bam.baseName
    def file_tag_new = file_name+'_BQSRecalibrated'
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
    def file_tag_new = file_name+'_BQSRecalibrated'
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


