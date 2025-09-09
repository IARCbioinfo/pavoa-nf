/*
========================================================================================
    MUTECT2 MODULE
========================================================================================
*/


process make_bed {

    input:
        path indexes
        val nsplit

    output:
        path '*_regions.bed'

    shell:
        """
        cat *fai | awk '{print \$1"	"0"	"\$2 }' | grep -v -P "alt|random|Un|chrEBV|HLA" > temp.bed
        grep -v '^track' temp.bed | sort -T \$PWD -k1,1 -k2,2n | bedtools merge -i stdin | awk '{print \$1" "\$2" "\$3}' | cut_into_small_beds.r $nsplit
        """

    stub:
        """
        touch chr1_regions.bed
        touch chr2_regions.bed
        touch chr3_regions.bed
        touch chrX_regions.bed
        """
}



process mutect {
    
    memory params.mem+'GB' 
    cpus params.cpu

    input:
        tuple val(file_tag), path(bamT), path(baiT), path(bamN), path(baiN), path(bed)
        path(ref)
        path(indexes)
        tuple path(PON), path(PON_tbi)
        tuple path(known_sites), path(known_sites_tbi)

    output:
        tuple val(file_tag), path("${printed_tag}_*.vcf"), path("${printed_tag}*stats*"), emit: calls
        tuple val(file_tag), path("${printed_tag}_f1r2.tar.gz"), emit : f1r2

    shell:
        def bed_tag = bed.baseName.replaceAll("[^a-zA-Z0-9 _-]+","")
        def printed_tag = "${file_tag}_" + bed_tag
        def input_t = "-I " + bamT.join(" -I ")
        def input_n = (bamN.baseName == 'None') ? "" : "-I ${bamN} -normal \$normal_name"
        def KS_option = known_sites ? "--germline-resource ${known_sites.get(0)}" : ""
        def PON_option = PON ? "--panel-of-normals ${PON.get(0)}" : ""
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

    publishDir "${params.output_folder}/mutect2/", mode: 'copy', saveAs: {filename ->
        if (filename.indexOf(".stats") > 0) "stats/$filename"
        else if (filename.indexOf(".vcf") > 0) "intermediate_calls/raw_calls/$filename"
    }
    
    input:
        tuple val(file_tag), path(vcf_files), path(txt_files)

    output:
        tuple val(file_tag), path("${file_tag}_calls.vcf"), path("${file_tag}_calls.vcf.stats"), emit : calls

    shell:
        input_stats = "-stats " + txt_files.join(" -stats ")
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
        gatk MergeMutectStats $input_stats -O ${file_tag}_calls.vcf.stats
        """

    stub:
        """
        touch ${file_tag}_calls.vcf ${file_tag}_calls.vcf.stats
        """
}

process ReadOrientationLearn {
            
    publishDir "${params.output_folder}/mutect2/stats/", mode: 'copy'

    input:
        tuple val(file_tag), path(f1r2_files)

    output:
        tuple val(file_tag), path("*model.tar.gz"), emit : ROmodel

    shell:
        input_f1r2 = "-I " + f1r2_files.join(" -I ")
        """
        gatk LearnReadOrientationModel ${input_f1r2} -O ${file_tag}_read-orientation-model.tar.gz
        """
    
    stub:
        """
        touch ${file_tag}_read-orientation-model.tar.gz
        """
}


process ContaminationEstimation {
   	
	cpus '16'
    memory '164 GB'
	    
    publishDir "${params.output_folder}/mutect/contamination/", mode: 'copy'

	input:
        tuple val(file_tag), path(bamT), path(baiT), path(bamN), path(baiN)
        path(ref)
        path(indexes)
        tuple path(snp_contam), path(snp_contam_tbi)

	output:
	    tuple val(file_tag), path("*contamination.table")
	
    shell:
	    """
        gatk --java-options "-Xmx${params.mem}G" GetPileupSummaries -R ${ref} -I $bamT -V $snp_contam -L $snp_contam -O ${bamT.baseName}_pileups.table
        gatk --java-options "-Xmx${params.mem}G" GetPileupSummaries -R ${ref} -I $bamN -V $snp_contam -L $snp_contam -O ${bamN.baseName}_pileups.table

	    gatk --java-options "-Xmx${params.mem}G" CalculateContamination -I ${bamT.baseName}_pileups.table ${bamN.baseName}_pileups.table -O ${bamT.baseName}_contamination.table
	    """
    
    stub:
        """
        touch ${pileupT.baseName}_contamination.table
        """
}


process FilterMuTectOutputs {

    publishDir "${params.output_folder}/mutect2/intermediate_calls/filtered", mode: 'copy'

    input:
        tuple val(tumor_normal_tag), path(vcf), path(stats), path(ROmodel), path(contam_tables)
        tuple path(fasta_ref), path(fasta_ref_fai), path(fasta_ref_gzi), path(fasta_ref_dict)

    output:
        tuple val(tumor_normal_tag), path("*filtered.vcf*")

    shell:
        RO = (ROmodel.baseName=="NO_ROmodel") ? "": "--ob-priors " + ROmodel.join(" --ob-priors ")
        contam = (contam_tables.baseName == "NO_contam") ? "" : "--contamination-table " + contam_tables.join(" --contamination-table ")
        """
        gatk FilterMutectCalls -R $fasta_ref -V $vcf $contam $RO -O ${tumor_normal_tag}_filtered.vcf
        """
    
    stub:
        """
        touch ${tumor_normal_tag}_filtered.vcf
        """
}


process FilterMuTectOutputsOnPass {

    publishDir "$params.output_folder/mutect/PASS/", mode: 'copy'

    input:
        tuple val(tumor_normal_tag), path(vcf_filtered)

    output:
        path("*_PASS.vcf*")

    shell:
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
    known_sites
    snp_contam

    main:

    // Panel Of Normal
    PON = params.PON ? (tuple file(params.PON), file(params.PON +'_TBI')) : (tuple file("NO_FILE"), file("NO_FILE_TBI"))

    //mutect2
    regions = make_bed(indexes,params.nsplit) | flatten
    mutect(pairs.combine(regions), ref, indexes, PON, known_sites)
    mutectOutput = mergeMuTectOutputs(mutect.out.calls.groupTuple(size: params.nsplit))

    // Read orientation
    ReadOrientationLearn(mutect.out.f1r2.groupTuple(size: params.nsplit))
    mutectOutput = mutectOutput.join( ReadOrientationLearn.out )

    // Estimate contamination
    if(snp_contam){
        contamination = ContaminationEstimation(pairs,ref,indexes, snp_contam)
        mutectOutput = mutectOutput.join(contamination)
    } else {
        mutectOutput = mutectOutput.map{row -> [row[0], row[1], row[2] , row[3], file("NO_contam") ]}
    }
    
    // filter
    FilterMuTectOutputs(mutectOutput, ref) | FilterMuTectOutputsOnPass

    emit:
    vcfs = FilterMuTectOutputs.out.calls

}