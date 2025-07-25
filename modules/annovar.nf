/*
========================================================================================
    ANNOTATION MODULE 
========================================================================================
*/

process annovar_annot {

    tag "${file_tag}"
    label 'annotation'

    memory = params.mem+'.GB'
    cpus params.cpu

    publishDir params.output_folder, mode: 'copy', pattern: '{*multianno*}'

    input:
        tuple val(file_tag), path(vcf)
        path annovar
        path annovarDB
        path annovarDBlist

    output:
        tuple val(file_tag), path("*avinput"), path("*multianno.txt"), path("*multianno.vcf"), emit: annotated

    shell:
        """
        annovar_annot.r -i ${vcf} -t ${params.cpu} -p "${params.pass}" \\
                -l ${annovarDBlist} -a ${annovarDB} -b ${annovar}
        for file in *multianno*; do
            mv \$file \${file/.vcf.gz/}
        done
        """

    stub:
        """
        touch ${file_tag}.tsv ${file_tag}_avinput ${file_tag}.multianno.txt ${file_tag}.multianno.vcf
        """

}


process gama_annot {

    tag "${file_tag}"
    label 'annotation'
    
    publishDir params.output_folder, mode: 'copy'

    memory = params.mem+'.GB'
    cpus params.cpu

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
        """
        touch ${file_tag}.1.tsv
        """   
}

process filter_vcf {

    tag "${file_tag}"
    label 'annotation'

    memory = params.mem+'.GB'
    cpus params.cpu

    publishDir "${params.output_folder}/filtered/${file_tag}/", mode: 'copy'

    input:
        tuple val(file_tag), path(tsv), path(vcf)

    output:
        tuple val(file_tag), path("${file_tag}/*vcf"), emit: filtered_vcf

    shell:
        """
        gama_filter.r $file_tag
        """

    stub:
        """
        mkdir -p ${file_tag}
        touch "${file_tag}/${file_tag}_snv.vcf"
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

    emit:
    annotated_tsv = filter_vcf.out.filtered_vcf


}