
include { DUPCALLER_INDEXES         } from './dupcaller'

workflow PREPARE_REFERENCES {

    take:
        ref_file
        bwa_index_dir
        output_dir

    main:
        ref = file(ref_file)
        assert (ref.extension == "fasta" || ref.extension == "fa") : "Error: reference file should be fasta format (.fasta or .fa)"

        ref_folder = ref.parent
        bwa_index_dir = (bwa_index_dir == null) ? ref_folder : bwa_index_dir
        ref_base = ref.name
        index_ch = Channel.empty()

        // Check if all BWA-mem2 index files exist
        required_extensions = ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac']
        index_files = required_extensions.collect { ext ->
            file(bwa_index_dir + "/${ref_base}.${ext}")
        }
        all_indexes_exist = index_files.every { it.exists() }

        // If all index files exist, use them; otherwise, create the index
        if (all_indexes_exist) {
            log.info(ref_base + ": bwa-mem2 indexes found")
            bwa_idx_ch = Channel.of( index_files.collect { file(it) })
            index_ch = index_ch.concat(bwa_idx_ch)
        } else {
            log.info(ref_base + ": bwa-mem2 indexes not found, creating them")
            bwa_idx_ch = CREATE_BWA_INDEXES(ref, output_dir).flatten()
            index_ch = index_ch.concat(bwa_idx_ch)
        }

        // Check if all dupcaller index files exist
        required_extensions = ['hp.h5', 'ref.h5', 'tn.h5']
        index_files = required_extensions.collect { ext ->
            file(bwa_index_dir + "/${ref_base}.${ext}")
        }
        all_indexes_exist = index_files.every { it.exists() }

        // If all index files exist, use them; otherwise, create the index
        if (all_indexes_exist) {
            log.info(ref_base + ": dupcaller indexes found")
            dup_idx_ch = Channel.of( index_files.collect { file(it) })
            index_ch = index_ch.concat(dup_idx_ch)
        } else {
            log.info(ref_base + ": dupcaller indexes not found, creating them")
            dup_idx_ch = DUPCALLER_INDEXES(ref, output_dir).flatten()
            index_ch = index_ch.concat(dup_idx_ch)
        }

        // Check if the .fai index exists
        sam_idx_path = file(ref_file + ".fai")
        if (sam_idx_path.exists()) {
            log.info(ref_base + ": samtools index found")
            sam_idx_ch = Channel.value(sam_idx_path)
            index_ch = index_ch.concat(sam_idx_ch)
        } else {
            log.info(ref_base + ": samtools index not found, creating it")
            sam_idx_ch = CREATE_SAMTOOLS_INDEXES(ref, output_dir)
            index_ch = index_ch.concat(sam_idx_ch)
        }

        // Check if the .dict index exists
        gatk_idx_path = file(ref_file.replaceAll(/\.(fa|fasta)$/, "") + ".dict")
        if (gatk_idx_path.exists()) {
            log.info(ref_base + ": GATK index found")
            gatk_idx_ch = Channel.value(gatk_idx_path)
            index_ch = index_ch.concat(gatk_idx_ch)
        } else {
            log.info(ref_base + ": GATK index not found, creating it")
            gatk_idx_ch = CREATE_GATK_INDEXES(ref, output_dir)
            index_ch = index_ch.concat(gatk_idx_ch)
        }

        // Check if the .alt index exists
        alt_idx_path = file(ref_file + ".alt")
        if(alt_idx_path.exists()) {
            log.info(ref_base + ": .alt index found")
            alt_idx_ch = Channel.value(alt_idx_path)
            index_ch = index_ch.concat(alt_idx_ch)
        }

    emit:
        ref = ref
        indexes = index_ch.flatten().collect()
}

workflow PREPARE_REFERENCES_DUPCALLER {

    take:
    ref_file
    output_dir

    main:

    ref = file(ref_file)
    assert (ref.extension == "fasta" || ref.extension == "fa") : "Error: reference file should be fasta format (.fasta or .fa)"

    ref_folder = ref.parent
    ref_base = ref.name
    index_ch = Channel.empty()

    // Check if all dupcaller index files exist
    required_extensions = ['hp.h5', 'ref.h5', 'tn.h5']
    index_files = required_extensions.collect { ext ->
        file(ref_folder + "/${ref_base}.${ext}")
    }
    all_indexes_exist = index_files.every { it.exists() }

    // If all index files exist, use them; otherwise, create the index
    if (all_indexes_exist) {
        log.info(ref_base + ": dupcaller indexes found")
        dup_idx_ch = Channel.of( index_files.collect { file(it) })
        index_ch = index_ch.concat(dup_idx_ch)
    } else {
        log.info(ref_base + ": dupcaller indexes not found, creating them")
        dup_idx_ch = DUPCALLER_INDEXES(ref, output_dir).flatten()
        index_ch = index_ch.concat(dup_idx_ch)
    }

    // Check if the .fai index exists
    sam_idx_path = file(ref_file + ".fai")
    if (sam_idx_path.exists()) {
        log.info(ref_base + ": samtools index found")
        sam_idx_ch = Channel.value(sam_idx_path)
        index_ch = index_ch.concat(sam_idx_ch)
    } else {
        log.info(ref_base + ": samtools index not found, creating it")
        sam_idx_ch = CREATE_SAMTOOLS_INDEXES(ref, output_dir)
        index_ch = index_ch.concat(sam_idx_ch)
    }

    emit:
    ref = ref
    indexes = index_ch.flatten().collect()
}


process CREATE_BWA_INDEXES {
    tag "${ref.baseName}"
    label 'alignment'

    publishDir "${output_dir}/index/", mode: 'copy'

    input:
    path(ref)
    val output_dir
    
    output:
    path("*.{0123,amb,ann,bwt.2bit.64,pac}")
    
    script:
    """
    # Create BWA-mem2 index
    bwa-mem2 index ${ref}
    """

    stub:
    """
    # Stub for BWA index creation
    touch ${ref.baseName}.0123 ${ref.baseName}.amb ${ref.baseName}.ann ${ref.baseName}.bwt.2bit.64 ${ref.baseName}.pac
    """
}


process CREATE_SAMTOOLS_INDEXES {
    tag "${ref.baseName}"
    label 'gatk'

    publishDir "${output_dir}/index", mode: 'copy'

    input:
    path(ref)
    val output_dir

    output:
    path("*.fai")

    script:
    """
    # Create .fai index with samtools
    samtools faidx ${ref}
    """

    stub:
    """
    touch ${ref}.fai
    """
}

process CREATE_GATK_INDEXES {
    tag "${ref.baseName}"
    label 'gatk'

    publishDir "${output_dir}/index", mode: 'copy'

    input:
    path(ref)
    val output_dir

    output:
    path("*.dict")

    script:
    """
    # Create .dict index with GATK
    gatk CreateSequenceDictionary \
        -R ${ref} \
        -O ${ref.baseName}.dict
    """

    stub:
    """
    touch ${ref.baseName}.dict
    """
}


process get_chromosomes {

    tag "$ref"

    input:
    path(ref)
    path(index)

    output:
    path("chrom_list.txt")

    script:
    """
    cut -f1 ${ref}.fai | sort | uniq | grep -v -P "alt|random|Un|chrEBV|HLA" > chrom_list.txt
    """

    stub:
    """
    cut -f1 ${ref}.fai | sort | uniq | grep -v -P "alt|random|Un|chrEBV|HLA" > chrom_list.txt
    """
}