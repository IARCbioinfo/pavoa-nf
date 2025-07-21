#!/usr/bin/env nextflow

process MERGE_BAM {
    tag "$file_tag"
    label 'alignment'
    
    cpus params.cpu
    memory params.mem+'G'

    if(!params.recalibration) {
        publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*.bam*"
    }

    input:
    tuple val(file_tag), path(bams), path(bais)

    output:
    tuple val(file_tag), path("${file_tag_new}.bam"), path("${file_tag_new}.bam.bai"), emit: merged_bam

    script:
    // Handle global parameters for mate tags
    bwa_opt = params.bwa_option_m ? '-M ' : ''
    samblaster_opt = params.bwa_option_m ? '-M ' : ''
    
    // Handle UMI and mode parameters for addMateTags
    addMateTags = ""
    addMateTags = params.umi ? "samblaster ${samblaster_opt} --addMateTags -a |" : "samblaster ${samblaster_opt} --addMateTags |"
    
    file_tag_new = file_tag
    if(params.trim) file_tag_new = file_tag_new + '_trimmed'
    if(params.alt) file_tag_new = file_tag_new + '_alt'
    
    
    if(bams instanceof List) {
        merge_threads = [params.cpu.intdiv(2) - 1, 1].max()
        sort_threads = [params.cpu.intdiv(2) - 1, 1].max()
        sort_mem = params.mem.div(2)
        bam_files = bams.join(" ")
        file_tag_new = file_tag_new + "_merged"
        
        """
        sambamba merge -t ${merge_threads} -l 0 /dev/stdout ${bam_files} | \\
        sambamba sort -t ${sort_threads} -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}.bam /dev/stdin
        """
    } else {
        """
        # Single BAM file - create symlinks with proper naming
        ln -sf ${bams[0]} ${file_tag_new}.bam
        ln -sf ${bais[0]} ${file_tag_new}.bam.bai
        """
    }

    stub:
    file_tag_new = file_tag
    if(params.trim) file_tag_new = file_tag_new + '_trimmed'
    if(params.alt) file_tag_new = file_tag_new + '_alt'
    if(bams.size() > 1) file_tag_new = file_tag_new + "_merged"
    """
    touch ${file_tag_new}.bam ${file_tag_new}.bam.bai
    """
}

