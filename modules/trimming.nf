/*
========================================================================================
    TRIMMING MODULE 
========================================================================================
*/


process trim_galore{

  tag "$file_tag"
  label 'trimgalore'

  cpus params.cpu

  publishDir "${params.output_folder}/QC/fastq/", mode: 'copy', pattern: '*.html,*.txt,*.zip'

  input:
  tuple val(file_tag), val(read_group), path(pair1), path(pair2)

  output:
  tuple val(file_tag), val(read_group), path("${file_tag}_trimmed_R1.fastq.gz"), path("${file_tag}_trimmed_R2.fastq.gz"), emit: trimmed_fastq
  path("*.{zip,txt}"), emit: reports
  path("*.html")

  script:
  """
  trim_galore --cores ${task.cpus} --${params.adapter} --quality ${params.quality} --length ${params.length} --fastqc --paired "$pair1" "$pair2" -o .
  mv ${file_tag}_1_val_1.fq.gz ${file_tag}_trimmed_R1.fastq.gz
  mv ${file_tag}_2_val_2.fq.gz ${file_tag}_trimmed_R2.fastq.gz
  mv ${file_tag}_1_val_1_fastqc.zip ${file_tag}_R1_fastqc.zip
  mv ${file_tag}_2_val_2_fastqc.zip ${file_tag}_R2_fastqc.zip
  mv ${file_tag}_1_val_1_fastqc.html ${file_tag}_R1_fastqc.html
  mv ${file_tag}_2_val_2_fastqc.html ${file_tag}_R2_fastqc.html
  """

  stub:
  """
  touch "${file_tag}_trimmed_R1.fastq.gz" "${file_tag}_trimmed_R2.fastq.gz"
  touch "${file_tag}_R1_fastqc.zip" "${file_tag}_R2_fastqc.zip"
  touch "${file_tag}_R1_fastqc.html" "${file_tag}_R2_fastqc.html"
  touch "${file_tag}_1_val_1.fastq.gz_trimming_report.txt" "${file_tag}_2_val_2.fastq.gz_trimming_report.txt"
  """
}


workflow TRIM_GALORE{

    take:
    samples

    main:
    // Read the input file and create a channel of tuples
    trim_galore(samples)

    emit:
    trimmed_fastq = trim_galore.out.trimmed_fastq
    reports = trim_galore.out.reports

}