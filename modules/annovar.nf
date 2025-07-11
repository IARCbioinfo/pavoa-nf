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

    output:
        tuple val(file_tag), path("*avinput"), path("*multianno.txt"), path("*multianno.vcf"), emit: annotated

    shell:
        """
        annovar_annot.r -i ${vcf} -t ${params.cpu} -p "${params.pass}" \\
                -l ${params.annovarDBlist} -a ${annovarDB} -b ${annovar}
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
        tuple val(file_tag), file("*1.tsv"), emit: annotated

    shell:
        """
        gama_annot.r -a ${annovarDB}
        """

    stub:
        """
        touch ${file_tag}.1.tsv
        """   
}


workflow ANNOTATION{

    take:
    vcfs

    main:
    def annovar = file( params.annovarBinPath )
    def annovarDB = file( params.annovarDBpath )

    annovar_annot(vcfs,annovar,annovarDB)
    gama_annot(annovar_annot.out.annotated,annovarDB)

    emit:
    annotated_tsv = gama_annot.out.annotated   

}