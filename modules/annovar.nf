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

    publishDir "${params.output_folder}/annotated/${file_tag}/", mode: 'copy'

    input:
        tuple val(file_tag), path(vcf)
        path annovar
        path annovarDB
        path annovarDBlist

    output:
        tuple val(file_tag), path("*_avinput"), path("*.txt"), path("*.vcf"), emit: annotated

    shell:
        def type = vcf.baseName.contains("indel") ? "_indel" : (vcf.baseName.contains("snv") ? "_snv" : "_call")
        """
        annovar_annot.r -i ${vcf} -t ${params.cpu} -p "${params.pass}" \\
                -l ${annovarDBlist} -a ${annovarDB} -b ${annovar}

        mv ${file_tag}*vcf ${file_tag}${type}.vcf
        mv ${file_tag}*avinput ${file_tag}${type}_avinput
        mv ${file_tag}*multianno*txt ${file_tag}${type}_multianno.txt
        """

    stub:
        def type = vcf.baseName.contains("indel") ? "_indel" : (vcf.baseName.contains("snv") ? "_snv" : "_call")
        """
        touch ${file_tag}${type}_avinput ${file_tag}${type}.txt ${file_tag}${type}.vcf
        """

}


process gama_annot {

    tag "${file_tag}"
    label 'annotation'
    
    publishDir "${params.output_folder}/annotated/${file_tag}/", mode: 'copy'

    memory 1+'GB'
    cpus 1

    input:
        tuple val(file_tag), path(avinput), path(tab), path(vcf)
        path annovarDB

    output:
        tuple val(file_tag), path("*1.tsv"), path(vcf), emit: annotated

    shell:
        """
        gama_annot.r -a ${annovarDB}
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

    memory = 1+'GB'
    cpus 1

    publishDir "${params.output_folder}/filtered/", mode: 'copy'

    input:
        tuple val(file_tag), path(tsv), path(vcf)

    output:
        tuple val(file_tag), path("${file_tag}/*vcf"), emit: filtered_vcf

    shell:
        //def new_tag = file_tag.replaceAll(/(\.[^.]+)?_multianno/, '')
        """
        gama_filter.r $file_tag ${params.cov_n_thresh} ${params.cov_t_thresh} ${params.min_vaf_t_thresh} ${params.max_vaf_t_thresh} ${params.cov_alt_t_thresh}
        #mv ${file_tag}/${file_tag}.vcf ${file_tag}/${file_tag}.vcf
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

    main:
    def annovar = file( params.annovarBinPath )
    def annovarDB = file( params.annovarDBpath )
    def annovarDBlist = file( params.annovarDBlist )

    annovar_annot(vcfs,annovar,annovarDB,annovarDBlist)
    gama_annot(annovar_annot.out.annotated,annovarDB)
    filter_vcf(gama_annot.out.annotated)
    filter_vcf.out.filtered_vcf.view()

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