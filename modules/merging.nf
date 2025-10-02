#!/usr/bin/env nextflow

process MERGE_BAM {
    tag "$file_tag"
    label 'alignment'
    
    cpus params.cpu
    memory params.mem+'G'

    publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*.bam*", saveAs: {filename ->
        if(!params.recalibration) filename
        else null
    }

    input:
    tuple val(file_tag), path(bams), path(bais)

    output:
    tuple val(file_tag), path("${file_tag}.bam"), path("${file_tag}.bam.bai"), emit: merged_bam

    script:
    if(bams instanceof List) {
        merge_threads = [params.cpu.intdiv(2) - 1, 1].max()
        sort_threads = [params.cpu.intdiv(2) - 1, 1].max()
        sort_mem = params.mem.div(2)
        bam_files = bams.join(" ")
        
        """
        sambamba merge -t ${merge_threads} -l 0 /dev/stdout ${bam_files} | \\
        sambamba sort -t ${sort_threads} -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag}.bam /dev/stdin
        """
    } else {
        """
        # Single BAM file - create symlinks with proper naming
        ln -sf ${bams[0]} ${file_tag}.bam
        ln -sf ${bais[0]} ${file_tag}.bam.bai
        """
    }

    stub:
    file_tag = file_tag
    """
    touch ${file_tag}.bam ${file_tag}.bam.bai
    """
}

