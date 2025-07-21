#!/usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

/*
========================================================================================
    ALIGNMENT MODULE - BWA-MEM ALIGNMENT PROCESSES
========================================================================================
*/

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process fastq_alignment {
    tag "$file_tag"
    label 'alignment'
    
    cpus params.cpu
    memory params.mem+'GB'    

    if(!params.recalibration){ 
        publishDir "${params.output_folder}/BAM/", mode: 'copy'	
    }

    input:
    tuple val(file_tag), val(read_group), path(pair1), path(pair2)
    path(ref)
    path(indexes)

    output:
    tuple val(file_tag), val(read_group), path("${file_tag_new}*.bam"), path("${file_tag_new}*.bai"), emit: bamfiles

    script:
    // Check if alt reference exists
    alt_aware = indexes.find { it.name.endsWith('.alt') }
    postaltjs_path = params.postaltjs ? params.postaltjs : "/opt/bwa-postalt.js"
    //postaltjs_path = "/opt/bwa-postalt.js"
    postalt = alt_aware ? 'k8 '+ postaltjs_path + ' ' + ref.name + '.alt' + ' | ' : ''
    
    // bwa user parameters
    bwa_opt = params.bwa_option_m ? '-M ' : ''

    // samblaster user parameters
    samblaster_opt = params.bwa_option_m ? '-M ' : ''

    // Handle UMI and mode parameters
    addMateTags = params.umi ? "samblaster ${samblaster_opt} --addMateTags -a |" : "samblaster ${samblaster_opt} --addMateTags |"
    
    // UMI -C option
    bwa_opt = params.umi ? bwa_opt + ' -C ' : bwa_opt
    
    // File naming
    pair2_input = (pair2.baseName == 'SINGLE') ? "" : pair2
    bwa_threads = [params.cpu.intdiv(2) - 1, 1].max()
    sort_threads = [params.cpu.intdiv(2) - 1, 1].max()
    sort_mem = [params.mem.intdiv(4), 1].max()
    
    file_tag_new = read_group == "" ? file_tag : "${file_tag}_${read_group}"
    if(alt_aware) file_tag_new = file_tag_new + '_alt'	
    
    RG = "\"@RG\\tID:${file_tag_new}\\tSM:${file_tag}\\t${params.pl}\""
    
    """
    set -o pipefail
    touch ${file_tag_new}.bam.bai
    bwa-mem2.sse42 mem ${bwa_opt} -t${bwa_threads} -R ${RG} ${ref} ${pair1} ${pair2_input} | \\
    ${postalt} ${addMateTags} \\
    sambamba view -S -f bam -l 0 /dev/stdin | \\
    sambamba sort -t ${sort_threads} -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}.bam /dev/stdin
    """

    stub:
    file_tag_new = read_group == "" ? file_tag : "${file_tag}_${read_group}"
    """
    touch ${file_tag_new}.bam ${file_tag_new}.bam.bai
    """
}

process bam_realignment {
    tag "$file_tag"
    label 'alignment'
    
    cpus params.cpu
    memory params.mem + 'G'
        
    if(!params.recalibration) {
        publishDir "${params.output_folder}/BAM/", mode: 'copy'
    }

    input:
    path file_tag
    tuple path(ref)
    tuple path(indexes)

    output:
    tuple val(file_tag), path("${file_tag_new}*.bam"), path("${file_tag_new}*.bai"), emit: realigned_bam

    script:
    // Check if alt reference exists
    alt_aware = indexes.find { it.name.endsWith('.alt') }
    postaltjs_path = params.postaltjs ? params.postaltjs : "/opt/conda/envs/alignment-nf/share/bwakit-0.7.15-1/bwa-postalt.js"
    postalt = alt_aware ? 'k8 '+ postaltjs_path + ' ' + ref.name + '.alt' + ' | ' : ''

    // bwa user parameters
    bwa_opt = params.bwa_option_m ? '-M ' : ''

    // samblaster user parameters
    samblaster_opt = params.bwa_option_m ? '-M ' : ''
    
    // UMI -C option
    bwa_opt = params.umi ? bwa_opt + ' -C ' : bwa_opt
    
    file_tag = file_tag.baseName
    file_tag_new = file_tag + '_realigned'
    if(alt_aware) file_tag_new = file_tag_new + '_alt'
    
    bwa_threads = [params.cpu.intdiv(2) - 1, 1].max()
    sort_threads = [params.cpu.intdiv(2) - 1, 1].max()
    sort_mem = params.mem.div(4)
    read_group = "@RG\\tID:${file_tag}\\tSM:${file_tag}\\t${params.pl}"
    
    """
    set -o pipefail
    
    samtools collate -uOn 128 ${infile} tmp_${file_tag} | \\
    samtools fastq - | \\
    bwa-mem2 mem ${bwa_opt} -t${bwa_threads} -R ${read_group} -p ${ref} - | \\
    ${postalt} samblaster ${samblaster_opt} --addMateTags --ignoreUnmated | \\
    sambamba view -S -f bam -l 0 /dev/stdin | \\
    sambamba sort -t ${sort_threads} -m ${sort_mem}G --tmpdir=${file_tag}_tmp -o ${file_tag_new}.bam /dev/stdin
    """

    stub:
    file_tag_new = file_tag + '_realigned'
    """
    touch ${file_tag_new}.bam ${file_tag_new}.bam.bai
    """
}


/*
========================================================================================
    WORKFLOWS
========================================================================================
*/

workflow FASTQ_ALIGNMENT {
    take:
    reads_ch        // channel: tuple(sample_id, nb_groups, read_group, fastq1, fastq2)
    reference_ch    // channel: tuple(ref files)
    indexes_ch      // channel: tuple(ref files)
    
    main:
    // Perform alignment
    fastq_alignment(reads_ch, reference_ch, indexes_ch)
    
    emit:
    bam_files = fastq_alignment.out.bamfiles
}

workflow BAM_REALIGNMENT {
    take:
    bam_files_ch    // channel: path(bam_file)
    reference_ch    // channel: tuple(ref files)
    indexes_ch      // channel: tuple(index files)
    
    main:
    bam_realignment(bam_files_ch, reference_ch, indexes_ch)

    emit:
    realigned_bam = bam_realignment.out.realigned_bam
}

workflow ALIGNMENT {
    take:
    input_ch        // channel: depends on mode (fastq or bam)
    reference_ch    // channel: tuple(ref files)
    indexes_ch      // channel: tuple(index files)

    main:
    //TODO: Determine mode based on input type
    //mode = input_ch.hasType(Tuple) ? 'fastq' : 'bam'
    if (mode == 'fastq') {
        FASTQ_ALIGNMENT(input_ch, reference_ch, indexes_ch)
        aligned_bam = FASTQ_ALIGNMENT.out.bam_files
    } else if (mode == 'bam') {
        BAM_REALIGNMENT(input_ch, reference_ch, indexes_ch)
        aligned_bam = BAM_REALIGNMENT.out.realigned_bam
    }
    
    emit:
    bam_files = aligned_bam
}
