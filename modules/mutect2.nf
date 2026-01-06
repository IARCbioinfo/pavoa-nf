/*
========================================================================================
    MUTECT2 MODULE
========================================================================================
*/


process make_bed {

    label 'mutect2'

    input:
    path bed
    path indexes
    val nsplit

    output:
    path '*_regions.bed'

    script:
    if(params.feature_file){
        """
        grep -v '^track' $bed | sort -T \$PWD -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print \$1" "\$2" "\$3}' | cut_into_small_beds.r $nsplit
        """
    } else {
        """
        cat *fai | awk '{print \$1"	"0"	"\$2 }' | grep -v -P "alt|random|Un|chrEBV|HLA" > temp.bed
        grep -v '^track' temp.bed | sort -T \$PWD -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print \$1" "\$2" "\$3}' | cut_into_small_beds.r $nsplit
        """
    }

}



process mutect {
    
    memory params.mem+'GB' 
    cpus params.cpu

    label 'mutect2'

    input:
        tuple val(file_tag), path(bamT), path(baiT), path(bamN), path(baiN), path(bed)
        path(ref)
        path(indexes)
        tuple path(PON), path(PON_tbi)
        path(known_sites)

    output:
        tuple val(file_tag), path("*_calls.vcf"), path("*_calls.vcf.stats"), emit: calls
        tuple val(file_tag), path("*_f1r2.tar.gz"), emit : f1r2

    script:
        def bed_tag = bed.baseName.replaceAll("[^a-zA-Z0-9 _-]+","")
        def printed_tag = "${file_tag}_" + bed_tag
        def input_t = "-I " + bamT.join(" -I ")
        def input_n = (bamN.baseName == 'None') ? "" : "-I ${bamN} -normal \$normal_name"
        def KS_vcf = known_sites.find { !it.name.endsWith(".tbi") }
        def KS_option = KS_vcf ? "--germline-resource " + KS_vcf.join(" --germline-resource ") : ""
        def PON_option = (PON.baseName == 'NO_FILE') ? "" : "--panel-of-normals ${PON.getAt(0)}"
        """
        normal_name=`samtools view -H $bamN | grep "^@RG" | head -1 | awk '{print \$NF}' | cut -c 4-`
        echo \$normal_name
        gatk Mutect2 --java-options "-Xmx${params.mem}G" -R $ref $KS_option $PON_option \
        $input_t $input_n -O ${printed_tag}_calls.vcf -L $bed $params.mutect_args --f1r2-tar-gz ${printed_tag}_f1r2.tar.gz
        """

    stub:
        def bed_tag = bed.baseName.replaceAll("[^a-zA-Z0-9 _-]+","")
        def printed_tag = "${file_tag}_" + bed_tag
        """
        touch ${printed_tag}_calls.vcf ${printed_tag}_calls.vcf.stats ${printed_tag}_calls_f1r2.tar.gz
        """
}

process mergeMuTectOutputs {

    label 'mutect2'

    publishDir "${params.output_folder}/mutect2/", mode: 'copy', saveAs: {filename ->
        if (filename.indexOf(".stats") > 0) "stats/$filename"
        else if (filename.indexOf(".vcf") > 0) "raw_calls/$filename"
    }
    
    input:
        tuple val(file_tag), path(vcf_files), path(txt_files)

    output:
        tuple val(file_tag), path("${file_tag}_calls.vcf"), path("${file_tag}_calls.vcf.stats"), emit : calls

    script:
        """
        # MERGE VCF FILES
        sed '/^#CHROM/q' `ls -1 *.vcf | head -1` > header.txt

        # Determine sort options based on the availability of version-sort in sort command
        sort_ops=\$(sort --help | grep -q 'version-sort' && echo "-k1,1V" || echo "-k1,1d")

        # Concatenate VCF contents and sort
        grep --no-filename -v '^#' *.vcf | LC_ALL=C sort -T \$PWD -t '	' \$sort_ops -k2,2n >> header.txt

        # Rename the merged file
        mv header.txt ${file_tag}_calls.vcf

        # MERGE STAT FILES
        
        # Build stats argument dynamically in bash
        input_stats=""
            for f in *stats; do
        input_stats="\$input_stats -stats \$f"
        done

        # Concatenate STAT files
        gatk MergeMutectStats \$input_stats -O ${file_tag}_calls.vcf.stats
        """

    stub:
        """
        touch ${file_tag}_calls.vcf ${file_tag}_calls.vcf.stats
        """
}

process ReadOrientationLearn {

    label 'mutect2'
            
    publishDir "${params.output_folder}/mutect2/stats/", mode: 'copy'

    input:
        tuple val(file_tag), path(f1r2_files)

    output:
        tuple val(file_tag), path("*model.tar.gz"), emit : ROmodel

    script:
        """
        # Build input_f1r2 argument dynamically in bash
        input_f1r2=""
        for f in *_f1r2.tar.gz; do
            input_f1r2="\$input_f1r2 -I \$f"
        done
        gatk LearnReadOrientationModel \${input_f1r2} -O ${file_tag}_read-orientation-model.tar.gz
        """
    
    stub:
        """
        touch ${file_tag}_read-orientation-model.tar.gz
        """
}


process ContaminationEstimation {
   	
	cpus '16'
    memory '164 GB'

    label 'mutect2'
	    
    publishDir "${params.output_folder}/mutect2/contamination/", mode: 'copy'

	input:
        tuple val(file_tag), path(bamT), path(baiT), path(bamN), path(baiN)
        path(ref)
        path(indexes)
        tuple path(snp_contam), path(snp_contam_tbi)

	output:
	    tuple val(file_tag), path("*contamination.table"), emit: contam_tables

    script:
	    """
        gatk --java-options "-Xmx${params.mem}G" GetPileupSummaries -R ${ref} -I $bamT -V $snp_contam -L $snp_contam -O ${bamT.baseName}_pileups.table
        gatk --java-options "-Xmx${params.mem}G" GetPileupSummaries -R ${ref} -I $bamN -V $snp_contam -L $snp_contam -O ${bamN.baseName}_pileups.table

	    gatk --java-options "-Xmx${params.mem}G" CalculateContamination -I ${bamT.baseName}_pileups.table ${bamN.baseName}_pileups.table -O ${bamT.baseName}_contamination.table
	    """
    
    stub:
        """
        touch ${bamN.baseName}_contamination.table
        """
}


process FilterMuTectOutputs {

    label 'mutect2'

    publishDir "${params.output_folder}/mutect2/calls/", mode: 'copy'

    input:
        tuple val(tumor_normal_tag), path(vcf), path(stats), path(ROmodel), path(contam_tables)
        path(ref)
        path(indexes)

    output:
        tuple val(tumor_normal_tag), path("*filtered.vcf{,.tbi}"), emit: calls

    script:
        RO = (ROmodel.baseName=="NO_ROmodel") ? "": "--ob-priors " + ROmodel.join(" --ob-priors ")
        contam = (contam_tables.baseName == "NO_contam") ? "" : "--contamination-table " + contam_tables.join(" --contamination-table ")
        """
        gatk FilterMutectCalls -R $ref -V $vcf $contam $RO -O ${tumor_normal_tag}_filtered.vcf
        """
    
    stub:
        """
        touch ${tumor_normal_tag}_filtered.vcf
        """
}


process FilterMuTectOutputsOnPass {

    label 'mutect2'

    publishDir "${params.output_folder}/mutect2/PASS/", mode: 'copy'

    input:
        tuple val(tumor_normal_tag), path(vcf_filtered)

    output:
        path("*_PASS.vcf*")

    script:
        vcf_name = vcf_filtered[0].baseName
        """
        bcftools view -f PASS -O z ${vcf_filtered[0]} > ${vcf_name}_PASS.vcf.gz
        bcftools index -t ${vcf_name}_PASS.vcf.gz
        """

    stub:
        vcf_name = vcf_filtered[0].baseName
        """
        touch ${vcf_name}_PASS.vcf.gz
        """
}

workflow MUTECT2_CALL{

    take:
    pairs // tuple (file_tag, bamT, baiT, bamN, baiN)
    ref
    indexes
    bed
    known_sites

    main:

    // Panel Of Normal
    PON = params.PON ? (tuple (file(params.PON), file(params.PON +'.tbi'))) : (tuple (file("NO_FILE"), file("NO_FILE_TBI")))

    //mutect2
    regions = make_bed(bed,indexes,params.nsplit) | flatten
    mutect(pairs.combine(regions), ref, indexes, PON, known_sites)
    mutectOutput = mergeMuTectOutputs(mutect.out.calls.groupTuple(size: params.nsplit))

    // Read orientation
    ReadOrientationLearn(mutect.out.f1r2.groupTuple(size: params.nsplit))
    mutectOutput = mutectOutput.join( ReadOrientationLearn.out )

    // Estimate contamination
    if(params.snp_contam){
        contam = tuple (file(params.snp_contam), file(params.snp_contam +'.tbi'))
        contamination = ContaminationEstimation(pairs,ref,indexes, contam)
        mutectOutput = mutectOutput.join(contamination)
    } else {
        mutectOutput = mutectOutput.map{row -> [row[0], row[1], row[2] , row[3], file("NO_contam") ]}
    }
    
    // filter
    FilterMuTectOutputs(mutectOutput, ref, indexes) //| FilterMuTectOutputsOnPass

    emit:
    // Collect all VCFs for annotation and add mutect2 tag
    vcfs = FilterMuTectOutputs.out.calls.map{ file_tag, vcf -> tuple(file_tag, vcf, "mutect2") }

}