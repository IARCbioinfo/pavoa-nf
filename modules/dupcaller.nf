/*
========================================================================================
    DUPCALLER MODULE 
========================================================================================
*/


process DUPCALLER_INDEXES {
    tag "${ref.baseName}"
    label 'dupcaller'

    publishDir "${output_dir}/index", mode: 'copy'

    input:
    path(ref)
    val output_dir

    output:
    path("*.h5")

    script:
    """
    DupCaller.py index -f ${ref}
    """

    stub:
    """
    touch "${ref.baseName}.hp.h5" "${ref.baseName}.ref.h5" "${ref.baseName}.tn.h5"
    """
}

process DUPCALLER_TRIM{

    tag "$file_tag"
    label 'dupcaller'

    input:
    tuple val(file_tag), val(read_group), path(pair1), path(pair2)

    output:
    tuple val(file_tag), val(read_group), path("${file_tag}_1.fastq.gz"), path("${file_tag}_2.fastq.gz"), emit: trimmed_fastq

    script:
    // Validate UMI tag parameter
    def umi = 'NNNNNNNN' // Default UMI tag if not specified
    umi = params.umi == true ? umi : params.umi.toUpperCase() // Use provided UMI tag if available

    """
    DupCaller.py trim -i "$pair1" -i2 "$pair2" -p "$umi" -o "$file_tag"
    # Compress output FASTQ files
    gzip "${file_tag}_1.fastq"
    gzip "${file_tag}_2.fastq"
    """

    stub:
    """
    touch "${file_tag}_1.fastq.gz" "${file_tag}_2.fastq.gz"
    """
}



process dupCallerCallByChromosome{

    tag "${file_tag}_${chromosome}"
    label 'dupcaller'

    memory params.mem +'GB'
    cpus params.cpu
    maxForks 15

    publishDir "${params.output_folder}/dupcaller/${chromosome}/", mode: 'copy'

    input:
        tuple val(file_tag), path(bamT), path(baiT), path(bamN), path(baiN)
        path(ref)
        path(indexes)
        path(known_sites)
        path(bed)
        each chromosome


    output:
        tuple val(file_tag), path("${file_tag}*/${file_tag}*_snv.vcf"), path("${file_tag}*/${file_tag}*_indel.vcf"), emit : vcfs
        path("${file_tag}*/*")

    script:
        def out_tag = "${file_tag}_${chromosome}"
        log.info "Processing ${out_tag}/${out_tag}_snv.vcf"
        def germline = known_sites.findAll { it.name.endsWith('.vcf.gz') }.collect { "-g ${it}" }.join(' ')
        def noise_mask = (bed.baseName=="NO_BED") ? "" : "-m " + bed.join(" -m ")
        """
        DupCaller.py call -tt 30 -b $bamT -n $bamN  -f $ref -o ${out_tag} -p $task.cpus -r $chromosome $germline $noise_mask
        DupCaller.py estimate -i ${out_tag} -f $ref -r $chromosome
        """

    stub:
        def out_tag = "${file_tag}_${chromosome}"
        """
        mkdir ${out_tag}
        cd ${out_tag}
        touch ${out_tag}_snv.vcf
        touch ${out_tag}_indel.vcf
        touch ${out_tag}.txt
        touch ${out_tag}.png
        """

}


process mergeResults{

    tag "${file_tag}"
    label 'dupcaller'

    memory params.mem+'GB'
    cpus params.cpu

    publishDir "${params.output_folder}/dupcaller/", mode: 'copy'

    input:
        tuple val(file_tag), path(snvs), path(indels)
        path(ref)
        path(indexes)

    output:
        tuple val(file_tag), path("${file_tag}/${file_tag}_calls.vcf"), emit: vcfs
        path("${file_tag}*/*")

    script:
        """
        # MERGE VCF FILES
        sed '/^#CHROM/q' `ls -1 *.vcf | head -1` > header.txt

        # Determine sort options based on the availability of version-sort in sort command
        sort_ops=\$(sort --help | grep -q 'version-sort' && echo "-k1,1V" || echo "-k1,1d")

        # Concatenate VCF contents and sort
        grep --no-filename -v '^#' *.vcf | LC_ALL=C sort -T \$PWD -t '	' \$sort_ops -k2,2n >> header.txt

        # Rename the merged file
        mkdir ${file_tag}
        mv header.txt ${file_tag}/${file_tag}_calls.vcf

        # Dupcaller Estimate
        DupCaller.py estimate -i ${file_tag} -f $ref
        """

    stub:
        """
        touch ${file_tag}_calls.vcf
        """
}


process dupCallerCall{

    tag "${file_tag}"
    label 'dupcaller'

    memory params.mem_dupcaller +'GB'
    cpus params.cpu_dupcaller

    publishDir "${params.output_folder}/dupcaller/", mode: 'copy'

    input:
        tuple val(file_tag), path(bamT), path(baiT), path(bamN), path(baiN)
        path(ref)
        path(indexes)
        path(known_sites)
        path(bed)
        val(mplDir)

    output:
        tuple val(file_tag), path("${file_tag}*/${file_tag}_snv.vcf"), 
                             path("${file_tag}*/${file_tag}_indel.vcf"),
                             path("${file_tag}*/*trinuc*"), emit: calls
        tuple val(file_tag), path("${file_tag}*/${file_tag}_snv.vcf"), emit : snvs
        tuple val(file_tag), path("${file_tag}*/${file_tag}_indel.vcf"), emit : indels
        tuple val(file_tag), path("${file_tag}*/*trinuc*"), emit: trinuc
        path("${file_tag}*/*")

    script:
        def germline = known_sites.findAll { it.name.endsWith('.vcf.gz') }.collect { "-g ${it}" }.join(' ')
        def noise_mask = (bed.baseName=="NO_BED") ? "" : "-m " + bed.join(" -m ")
        def chromosome = (params.chromosome) ? "-r ${params.chromosome}" : ""
        """
        export MPLCONFIGDIR="$mplDir"
        DupCaller.py call -tt 30 -b $bamT -n $bamN -f $ref -o ${file_tag} -p $task.cpus $germline $noise_mask ${chromosome}
        DupCaller.py estimate -i ${file_tag} -f $ref ${chromosome}
        """

    stub:
        """
        mkdir ${file_tag}
        cd ${file_tag}
        touch ${file_tag}_snv.vcf
        touch ${file_tag}_indel.vcf
        touch ${file_tag}.txt
        touch ${file_tag}.png
        touch ${file_tag}_trinuc_by_duplex_group.txt
        """

}


process fixDupcallerOutput{

    tag "${file_tag}"
    label 'dupcaller'

    memory params.mem+'GB'
    cpus params.cpu

    input:
        tuple val(file_tag), path(calls)

    output:
        tuple val(file_tag), path("${file_tag}_*.fixed.vcf.gz"), emit: vcfs

    script:
        def out = calls.baseName + ".fixed.vcf"
        """
        # add GT tag for annovar
        first_format_num=\$(grep -n -m 1 '##FORMAT' "$calls" | cut -d : -f 1)
        sed "\${first_format_num}a##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">" "$calls" > "$out"
        sed -E 's/(AC:RC:DP[[:space:]]+)([^:]+:[^:]+:[^:]+)[[:space:]]+([^:]+:[^:]+:[^:]+)/GT:\\10\\/1:\\2\\t0\\/0:\\3/' "$calls" > "$out"
        gzip -f $out
        """

    stub:
        def out = calls.baseName + ".fixed.vcf"
        """
        touch "${out}.gz"
        """
}

process DUPCALLER_ESTIMATE{

    tag "${file_tag}"
    label 'dupcaller'

    memory 1+'GB'
    cpus 4

    publishDir "${params.output_folder}/filtered/", mode: 'copy'

    input:
        tuple val(file_tag), path(snvs), path(indels), path(trinuc)
        path(ref)
        path(indexes)

    output:
        path("${file_tag}/*")

    script:
        def chromosome = (params.chromosome) ? "-r ${params.chromosome}" : ""
        """
        mkdir -p "${file_tag}"
        mv ${file_tag}_* ${file_tag}/.
        DupCaller.py estimate -i ${file_tag} -f $ref ${chromosome}
        """

    stub:
        """
        mkdir -p "${file_tag}"
        touch "${file_tag}/${file_tag}.txt"
        touch "${file_tag}/${file_tag}.png"
        """
}


workflow DUPCALLER_CALL{

    take:
    pairs
    ref
    indexes
    known_sites
    noise_mask

    main:

    // Create MPLCONFIGDIR for matplotlib
    def mplDir = "${workflow.workDir}/mplconfig"
    new File(mplDir as String).mkdirs()

    dupCallerCall(pairs, ref, indexes, known_sites, noise_mask, mplDir)

    DUPCALLER_ESTIMATE(dupCallerCall.out.calls,ref,indexes)
    snvindels=dupCallerCall.out.snvs.mix(dupCallerCall.out.indels)
    fixDupcallerOutput(snvindels)
    
    emit:
    vcfs = fixDupcallerOutput.out.vcfs
    trinuc = dupCallerCall.out.trinuc

}
