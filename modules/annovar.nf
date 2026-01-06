/*
========================================================================================
    ANNOTATION MODULE 
========================================================================================
*/

process annovar_annot {

    tag "${file_tag}"
    label 'annotation'

    memory params.mem+'.GB'
    cpus params.cpu

    publishDir "${params.output_folder}/${caller_tag}/calls/annovar/", mode: 'copy', pattern: "*multianno*"

    input:
        tuple val(file_tag), path(vcf), val(caller_tag)
        path annovar
        path annovarDB
        path annovarDBlist

    output:
        tuple val(file_tag), path("*multianno.txt"), path("${file_tag}/*vcf*"), val(caller_tag), emit: annotated
        path("*multianno*"), emit: multianno

    script:
        def type = vcf.baseName.contains("indel") ? "_indel" : (vcf.baseName.contains("snv") ? "_snv" : "_call")
        def db = annovarDB.baseName.replaceAll(/db/, '')
        """
        mkdir $file_tag
        cp -L ${vcf} ${file_tag}/${vcf.getName()}
        annovar_annot.r -i ${vcf} -t ${params.cpu} -p "${params.pass}" \\
                -l ${annovarDBlist} -a ${annovarDB} -b ${annovar}
        mv ${file_tag}*multianno.vcf ${file_tag}${type}.${db}_multianno.vcf
        mv ${file_tag}*multianno.txt ${file_tag}${type}.${db}_multianno.txt
        """

    stub:
        def type = vcf.baseName.contains("indel") ? "_indel" : (vcf.baseName.contains("snv") ? "_snv" : "_call")
        def db = annovarDB.baseName.replaceAll(/db/, '')
        """
        touch ${file_tag}${type}.${db}_multianno.vcf ${file_tag}${type}.${db}_multianno.txt
        """

}


process gama_annot {

    tag "${file_tag}"
    label 'annotation'
    
    publishDir "${params.output_folder}/${caller_tag}/calls/annovar/", mode: 'copy', pattern: "*1.tsv*"

    memory '1.GB'
    cpus 1

    input:
        tuple val(file_tag), path(multianno), path(vcf), val(caller_tag)
        path annovarDB
        path ref

    output:
        tuple val(file_tag), path("*1.tsv"), path(vcf), val(caller_tag), emit: annotated

    script:
        """
        gama_annot.r -a ${annovarDB} -r ${ref}
        """

    stub:
        def type = vcf.baseName.contains("indel") ? "_indel" : (vcf.baseName.contains("snv") ? "_snv" : "_call")
        """
        touch ${file_tag}${type}.1.tsv
        """   
}

process filter_vcf {

    tag "${file_tag}"
    label 'annotation'

    memory '1.GB'
    cpus 1

    publishDir "${params.output_folder}/${caller_tag}/calls/", mode: 'copy'

    input:
        tuple val(file_tag), path(tsv), path(vcf), val(caller_tag)

    output:
        tuple val(file_tag), path("filtered.1/*vcf"), val(caller_tag), emit: filtered_vcf

    script:
        """
        gama_filter.r $file_tag ${params.cov_n_thresh} ${params.cov_t_thresh} \
         ${params.min_vaf_t_thresh} ${params.max_vaf_t_thresh} \
         ${params.min_cov_alt_t_thresh} ${params.max_cov_alt_t_thresh} 
        """

    stub:
        def type = vcf.baseName.contains("indel") ? "_indel" : (vcf.baseName.contains("snv") ? "_snv" : "_call")
        """
        mkdir -p ${file_tag}
        touch "${file_tag}/${file_tag}${type}.vcf"
        """
}

workflow ANNOTATION{

    take:
    vcfs
    ref

    main:
    def annovar = file( params.annovarBinPath )
    def annovarDB = file( params.annovarDBpath )
    def annovarDBlist = file( params.annovarDBlist )

    annovar_annot(vcfs,annovar,annovarDB,annovarDBlist)
    gama_annot(annovar_annot.out.annotated,annovarDB, ref)
    filter_vcf(gama_annot.out.annotated)

    emit:
    annotated_tsv = filter_vcf.out.filtered_vcf

}


workflow  GAMA_FILTER{
    take:
    file_pairs

    main:

    // Apply GAMA filtering
    filter_vcf(file_pairs)

    emit:
    filtered_vcf = filter_vcf.out.filtered_vcf
}