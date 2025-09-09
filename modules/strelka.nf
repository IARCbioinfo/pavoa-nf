/*
========================================================================================
    STRELKA MODULE
========================================================================================
*/


process strelka_somatic {

    cpus params.cpu
    memory params.mem+'GB' 
      
    publishDir "${params.output_folder}/strelka/", mode: 'copy', pattern: "*vcf*"
    publishDir "${params.output_folder}/strelka/CallableRegions/", mode: 'copy', pattern: "*bed*"

    input:
        tuple val(file_tag), path(bamT), path(baiT), path(bamN), path(baiN)
        path(ref)
        path(indexes)
        path(bed)

     output:
        tuple val(file_tag), path("*.somatic.snvs.vcf.gz"), path("*.somatic.indels.vcf.gz"), emit: calls
        tuple val(file_tag), path("*.somatic.snvs.vcf.gz"), emit : snvs
        tuple val(file_tag), path("*.somatic.indels.vcf.gz"), emit : indels
        path("*callable.regions.bed.gz*"), emit : regionfiles
        path("*.somatic.*")

    script:
    def strelka = "configureStrelkaSomaticWorkflow.py"
    def config = params.strelka_config ? params.strelka_config : strelka + ".ini"
    def exome = params.exome ? "--exome" : ""
    def callRegions = (bed.baseName=="NO_BED") ? "" : "--callRegions " + bed

    """
    ${strelka} --tumorBam ${bamT} --normalBam ${bamN} --referenceFasta ${ref} \\
        --config ${config} ${exome} ${callRegions} \\
        --runDir strelkaAnalysis --outputCallableRegions
    
    ./strelkaAnalysis/runWorkflow.py -m local -j ${params.cpu} -g ${params.mem}

    mv strelkaAnalysis/results/variants/* .
    mv somatic.indels.vcf.gz !{file_tag}.indels.vcf.gz
    mv somatic.snvs.vcf.gz !{file_tag}.snvs.vcf.gz
    mv somatic.indels.vcf.gz.tbi !{file_tag}.indels.vcf.gz.tbi
    mv somatic.snvs.vcf.gz.tbi !{file_tag}.snvs.vcf.gz.tbi
    fixStrelkaOutput.sh *.vcf.gz

    mv strelkaAnalysis/results/regions/* .
    mv somatic.callable.regions.bed.gz !{file_tag}.callable.regions.bed.gz
    mv somatic.callable.regions.bed.gz.tbi !{file_tag}.callable.regions.bed.gz.tbi
    """
}


workflow STRELKA2_CALL{

    take:
    pairs
    ref
    indexes
    bed

    main:
    strelka_somatic(pairs, ref, indexes, bed)
    snvindels=strelka_somatic.out.snvs.mix(strelka_somatic.out.indels)
    
    emit:
    vcfs = snvindels

}
