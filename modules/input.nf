/*INPUT PREPARATION MODULE
========================================================================================
*/

/*
========================================================================================
    PROCESSES
========================================================================================
*/

// Process to validate input files exist and are readable
process VALIDATE_INPUT_FILES {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), val(read_group), path(pair1), path(pair2), val(normal)

    output:
    tuple val(sample_id), val(read_group), path(pair1), path(pair2), emit: validated_files

    script:
    def pair2_check = pair2.name != 'SINGLE' ? """
        # Check pair2 file
        if [[ ! -f "${pair2}" || ! -r "${pair2}" ]]; then
            echo "ERROR: Read file ${pair2} not found or not readable"
            exit 1
        fi
    """ : ""
    
    """
    # Check pair1 file
    if [[ ! -f "${pair1}" || ! -r "${pair1}" ]]; then
        echo "ERROR: Read file ${pair1} not found or not readable"
        exit 1
    fi
    
    ${pair2_check}
    
    echo "Input validation completed for sample ${sample_id}"
    """

    stub:
    """
    echo "Input validation completed for sample ${sample_id} (stub)"
    """
}

// Process to validate BAM files
process VALIDATE_BAM_FILES {
    tag "${bam.baseName}"
    
    input:
    path bam

    output:
    path bam, emit: validated_bam

    script:
    """
    # Check BAM file exists and is readable
    if [[ ! -f "${bam}" || ! -r "${bam}" ]]; then
        echo "ERROR: BAM file ${bam} not found or not readable"
        exit 1
    fi
    
    # Check if file is a valid BAM file
    if ! samtools quickcheck "${bam}"; then
        echo "ERROR: ${bam} is not a valid BAM file"
        exit 1
    fi
    
    echo "BAM validation completed for ${bam}"
    """

    stub:
    """
    echo "BAM validation completed for ${bam} (stub)"
    """
}

// Process to check input file format and headers
process VALIDATE_INPUT_TSV {
    tag "${input_file.baseName}"
    
    input:
    path input_file

    output:
    path input_file, emit: validated_tsv

    script:
    """
    # Check if file exists
    if [[ ! -f "${input_file}" ]]; then
        echo "ERROR: Input file ${input_file} not found"
        exit 1
    fi
    
    # Check if file has content
    if [[ ! -s "${input_file}" ]]; then
        echo "ERROR: Input file ${input_file} is empty"
        exit 1
    fi
    
    # Check for required columns (assuming header exists)
    if ! head -1 "${input_file}" | grep -q "SM"; then
        echo "ERROR: Missing required column 'SM' in ${input_file}"
        exit 1
    fi
    
    if ! head -1 "${input_file}" | grep -q "RG"; then
        echo "ERROR: Missing required column 'RG' in ${input_file}"
        exit 1
    fi
    
    if ! head -1 "${input_file}" | grep -q "pair1"; then
        echo "ERROR: Missing required column 'pair1' in ${input_file}"
        exit 1
    fi
    
    # Count number of data rows (excluding header)
    data_rows=\$(tail -n +2 "${input_file}" | wc -l)
    if [[ \$data_rows -eq 0 ]]; then
        echo "ERROR: No data rows found in ${input_file}"
        exit 1
    fi
    
    echo "Input TSV validation completed for ${input_file}"
    echo "Found \$data_rows samples to process"
    """

    stub:
    """
    echo "Input TSV validation completed for ${input_file} (stub)"
    """
}

/*
========================================================================================
    WORKFLOWS
========================================================================================
*/

workflow PREPARE_INPUT_FROM_TSV {
    take:
    input_file_path     // path: Input TSV file with sample information
    
    main:
    
    // Validate input TSV file first
    VALIDATE_INPUT_TSV(file(input_file_path))
    
    // Parse the input file and create sample channel
    readPairs = Channel.fromPath(input_file_path)
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row ->
            // Create tuple with sample info
            tuple(
                row.SM,                                       // Sample ID
                row.RG,                                       // Read Group ID
                file(row.pair1),                              // First pair file
                row.pair2 ? file(row.pair2) : file("SINGLE"), // Second pair file or SINGLE
                row.sample ? file(row.sample) : "NO_NORMAL"            // normal sample (optional, can be empty)
            )
        }
    
    // Validate each input file
    VALIDATE_INPUT_FILES(readPairs)
    
    emit:
    reads = VALIDATE_INPUT_FILES.out.validated_files
}



workflow PREPARE_INPUT_FROM_FOLDER {
    take:
    input_folder_path  // path: Folder containing FASTQ files
    fastq_ext          // val: FASTQ file extension
    suffix1            // val: Suffix for first pair
    suffix2            // val: Suffix for second pair
    
    main:
    
    // Check if folder contains FASTQ files
    fastq_pattern = "${input_folder_path}/*{${suffix1},${suffix2}}.${fastq_ext}"
    
    readPairs = Channel.fromFilePairs(fastq_pattern)
        .ifEmpty { error "No FASTQ files found in ${input_folder_path} with pattern *{${suffix1},${suffix2}}.${fastq_ext}" }
        .map { row -> 
            tuple(
                row[0],      // Sample ID (file_tag)
                "NO_RG",     // Read Group (read_group)
                row[1][0],   // Pair1 file
                row[1][1]    // Pair2 file
            )
        }
    
    // Validate input files
    VALIDATE_INPUT_FILES(readPairs)
    
    emit:
    reads = VALIDATE_INPUT_FILES.out.validated_files
}

workflow PREPARE_BAM_INPUT {
    take:
    input_folder_path   // path: Folder containing BAM files
    
    main:
    
    // Find BAM files in input folder
    bam_files = Channel.fromPath("${input_folder_path}/*.bam")
        .ifEmpty { error "No BAM files found in ${input_folder_path}" }
    
    // Validate BAM files
    VALIDATE_BAM_FILES(bam_files)
    
    emit:
    bams = VALIDATE_BAM_FILES.out.validated_bam
}

workflow PREPARE_CALLING_INPUT {

    take:
    input_file_path
    alignments // tuple(sample, bam, bai)

    main:

    //alignments.view()
    alignments_ch = alignments.flatten().collate(3).map { it -> tuple(it[0], it[1], it[2]) }
    //alignments_ch.view()

    // Channel of (tumor_sample, normal_sample)
    sample_pairs = Channel
        .fromPath(input_file_path)
        .splitCsv(header: true, sep: '\t', strip: true)
        .filter { row -> row.normal } // Only keep samples with matched normals
        .map { row -> tuple(row.SM, row.normal) }
        .unique()
    //sample_pairs.view()

    // First join: get tumor sample alignments
    tumor_aligned = sample_pairs
        .join(alignments_ch, by: 0) // match on tumor sample
        .filter { tumor_sample, normal_sample, _tumor_bam, _tumor_bai -> 
            tumor_sample != normal_sample // for stupid users that might put the same sample as both tumor and normal !
        }
        .map { tumor_sample, normal_sample, tumor_bam, tumor_bai ->
            tuple(normal_sample, tumor_sample, tumor_bam, tumor_bai)
        }
    //tumor_aligned.view()


    // Second join: get normal sample alignments
    full_pairs = tumor_aligned
        .combine(alignments_ch, by: 0) // match on normal sample
        .map { _normal_sample, tumor_sample, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(tumor_sample, tumor_bam, tumor_bai, normal_bam, normal_bai)
        }

    //full_pairs.view()

    emit:
    pairs = full_pairs
}
