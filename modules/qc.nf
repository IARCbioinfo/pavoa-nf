/*
========================================================================================
    QUALITY CONTROL MODULE
========================================================================================
*/

process MULTIQC {
    
    tag "multiqc"
    label 'multiqc'

    cpus params.cpu
    memory "4.GB"
    
    publishDir "${params.output_folder}/QC/", 
               mode: 'copy', 
               pattern: "multiqc_report.html"
    
    publishDir "${params.output_folder}/QC/", 
               mode: 'copy', 
               pattern: "multiqc_report_data/"
    
    input:
    path(reports)
    path(multiqc_config)
    
    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data/"), emit: data
    
    script:
    config = (multiqc_config.name=="NO_FILE") ? "" : "--config ${multiqc_config}"
    """
    set -euo pipefail
    
    # Create MultiQC report
    multiqc --force --filename multiqc_report.html ${config} .
    """

    stub:
    """
    mkdir -p multiqc_report_data && touch multiqc_report.html
    """
}

process FLAGSTAT {

    tag "$file_tag"
    label 'gatk'

    cpus params.cpu
    memory '1 G'

    publishDir "${params.output_folder}/QC/BAM/flagstat/", mode: 'copy'
    
    input:
    tuple val(file_tag), path(bam), path(bai)
    
    output:
    path("*.stats.txt"), emit: report
    
    script:
    """
    set -euo pipefail
    
    # Create flagstat report
    sambamba flagstat -t ${task.cpus} $bam > ${file_tag}.stats.txt
    """

    stub:

    """
    touch ${file_tag}.stats.txt
    """
}

process QUALIMAP {
    
    tag "$file_tag"
    label 'qualimap'

    cpus params.cpu
    memory params.mem + 'G'

    publishDir "${params.output_folder}/QC/BAM/qualimap/", mode: 'copy', pattern: "$file_tag/*"
    
    input:
    tuple val(file_tag), path(bam), path(bai)
    path(feature_file)
    
    output:
    path(file_tag), emit: report
    
    script:
    def feature = feature_file.name != 'NO_FILE' ? "--feature-file $feature_file" : ''
    """
    set -euo pipefail
    
    # Create qualimap report
    qualimap bamqc -nt $params.cpu $feature --skip-duplicated -bam $bam --java-mem-size=${params.mem}G -outdir $file_tag -outformat html
    """
    
    stub:
    """
    mkdir -p "${file_tag}" && touch "${file_tag}/report.html" "${file_tag}/report.txt"
    """
}

