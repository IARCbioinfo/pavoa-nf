/*
========================================================================================
    TRIMMING MODULE 
========================================================================================
*/


process trim_galore{

  tag "$file_tag"
  label 'trimgalore'

  cpus params.cpu

  publishDir "${params.output_folder}/QC/fastq/", mode: 'copy', pattern: ~/.*\.(html|txt|zip)$/

  input:
  tuple val(file_tag), val(read_group), path(pair1), path(pair2)

  output:
  tuple val(file_tag), val(read_group), path("*_trimmed_1.fq.gz"), path("*_trimmed_2.fq.gz"), emit: trimmed_fastq
  path("*.{zip,txt}"), emit: reports
  path("*.html")

  script:
  def basename1 = pair1.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
  def basename2 = pair2.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
  """
  trim_galore --cores ${task.cpus} --${params.adapter} --quality ${params.quality} --length ${params.length} --fastqc --paired "$pair1" "$pair2" -o .
  mv "${basename1}_val_1.fq.gz" "${file_tag}_${read_group}_trimmed_1.fq.gz"
  mv "${basename2}_val_2.fq.gz" "${file_tag}_${read_group}_trimmed_2.fq.gz"
  mv "${basename1}_trimming_report.txt" "${file_tag}_${read_group}_trimming_report.txt"
  mv "${basename2}_trimming_report.txt" "${file_tag}_${read_group}_trimming_report.txt"
  mv "${basename1}_val_1_fastqc.zip" "${file_tag}_${read_group}_1_fastqc.zip"
  mv "${basename2}_val_2_fastqc.zip" "${file_tag}_${read_group}_2_fastqc.zip"
  mv "${basename1}_val_1_fastqc.html" "${file_tag}_${read_group}_1_fastqc.html"
  mv "${basename2}_val_2_fastqc.html" "${file_tag}_${read_group}_2_fastqc.html"
  """

  stub:
  """
  touch "${file_tag}_${read_group}_trimmed_1.fq.gz" "${file_tag}_${read_group}_trimmed_2.fq.gz"
  touch "${file_tag}_${read_group}_1_fastqc.zip" "${file_tag}_${read_group}_2_fastqc.zip"
  touch "${file_tag}_${read_group}_1_fastqc.html" "${file_tag}_${read_group}_2_fastqc.html"
  touch "${file_tag}_${read_group}_trimming_report.txt" "${file_tag}_${read_group}_trimming_report.txt"
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