# pavoa-nf

![pavoa logo](assets/pavoa.png)

## Pipeline for variant detection: Preprocessing, Alignment, Variant calling, Optimization, Annotation

## Description

pavoa-nf is a comprehensive Nextflow pipeline for variant detection from whole genome or whole exome sequencing data. The pipeline performs preprocessing, alignment, variant calling, optimization, and annotation steps to produce high-quality variant calls with comprehensive annotations.

The pipeline integrates industry-standard tools and follows GATK best practices for variant detection, providing a streamlined workflow from raw sequencing reads to annotated variants.

## Dependencies

Containers are available with all the tools needed to run the pipeline (see nextflow.config and Usage section). Only Annovar require a local installation.

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.

2. External software:
   - [Trim-galore](https://github.com/FelixKrueger/TrimGalore) read trimming 
   - [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) fast alignment
   - [samblaster](https://github.com/GregoryFaust/samblaster) fast and flexible program for marking duplicates
   - [sambamba](https://github.com/lomereiter/sambamba) fast processing of NGS alignment 
   - [GATK4](https://software.broadinstitute.org/gatk/guide/quickstart) GATK tools : MarkDuplicates, BaseRecalibrator
   - [Dupcaller](https://github.com/AlexandrovLab/DupCaller) UMI trimming and Caller for UDseq
   - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) read quality assessment
   - [Qualimap](http://qualimap.conesalab.org/) alignment quality control
   - [MultiQC](https://multiqc.info/) quality control reports
   - [Annovar](http://annovar.openbioinformatics.org/en/latest/user-guide/download/) Variant annotation

3. Optional external software for advanced features:
   - the k8 javascript execution shell (e.g., available in the [bwakit](https://sourceforge.net/projects/bio-bwa/files/bwakit/) archive); must be in the PATH
   - javascript bwa-postalt.js and the additional fasta reference .alt file from [bwakit](https://github.com/lh3/bwa/tree/master/bwakit) must be in the same directory as the reference genome file (for alternative contig handling)

4. Reference files:
   - You can generate indexes for bwa-mem2, GATK and dupcaller, and store them with the reference fasta file. 
   - Or you can let pavoa-nf generate them for you. They will be available in the output folder.

5. VCF files :
   - Lists of indels and SNVs, vcf files and corresponding tabix indexes (.tbi)
   - Recommended: af-only-gnomad.hg38.vcf.gz, Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.


## Parameters

### Mandatory Parameters

#### Inputs

| Type | Description |
|------|-------------|
| --input_file | Input file (comma-separated) with 5 columns: SM (sample name), RG (read group ID), pair1 (first fastq of the pair), pair2 (second fastq of the pair), normal (SM of the normal sample for somatic calling) |

#### Reference

| Type  | Description |
|-------|-------------|
| --ref | genome reference with its index files for mapping (.fai, .ann, .amb, .pac, .bwt.2bit.64, .0123, .dict and optionnaly .alt), and for mapping (.hp.h5, .ref.h5 and .tn.h5); in the same directory) |

##### ALT contigs

Some regions in the human reference genome (e.g., GRCh38) have alternative haplotypes, represented as ALT contigs. These may cause ambiguous mappings when reads align equally well to both primary and ALT regions.

If a .alt index is provided, the pipeline automatically runs bwa-postalt.js to adjust mapping qualities and resolve ambiguities by linking ALT contigs to their primary counterparts.

##### Annovar DB

1- Install annovar executables (table_annovar.pl, ...). 
Add `--annovarDBpath /path/to/annovarbin/` to your pavoa-nf command.

2- Download the databases (eg. hg38db).
Add `--annovarDBpath /path/to/hg38db` to your pavoa-nf command.

3- To activate annotation step, provide a DBlist file (eg. [hg38_repeat_1.txt](assets/annovar/hg38db/hg38_repeat_1.txt) )
Add `--annovarDBlist hg38_repeat_1.txt` to your pavoa-nf command.


### Optional parameters

#### Outputs

| Name | Description |
|------|-------------|
| --output_folder | Output folder for results |

### Ressources

| Name | Default value | Description |
|------|---------------|-------------|
| --cpu | 8  | number of CPUs |
| --mem | 32 | memory (GB) |
| --cpu_bqsr | 2  | number of CPUs for GATK base quality score recalibration |
| --mem_bqsr | 10 | memory for GATK base quality score recalibration (GB) |


### Trimming options

| Name | Default value | Description |
|------|---------------|-------------|
| --trim | |enable adapter sequence trimming |
| --adapter | illumina | adapter type (illumina, nextera, etc.) |
| --length | 30 | Minimum read length after trimming |
| --quality | 30 | Minimum read quality after trimming |

### Processing Options

| Name | Default value | Description |
|------|---------------|-------------|
| --bqsr | | Enable base quality score recalibration |
| --umi | NNNNNNNN | Enable UMI-aware duplicate marking. You may pass a TAG or keep default if used as a flag |

### Quality Control

| Name | Description |
|------|-------------|
| --feature_file | Feature file (bed) for Qualimap |
| --multiqc_config | MultiQC configuration file |

### Known Sites and Masking

| Name | Description |
|------|-------------|
| --known_sites | VCF file, known variants to filter |
| --mask | BED file, regions to ignore during calling |

### Annotation

| Name | Description |
|------|-------------|
| --annovarDBlist | File with two columns : protocols and operations [see example](assets/annovar/hg38db/hg38_repeat_1.txt) |
| --annovarDBpath | Path to your annovarDB |
| --annovarBinPath | Path to table_annovar.pl |

## Usage

### Basic usage

To run the pipeline on a series of fastq files listed in input.txt and a fasta reference file `hg38.fasta`, one can type:

```bash
nextflow run iarcbioinfo/pavoa-nf -profile singularity --input_file input.txt --ref hg38.fasta
```

### Complete variant calling workflow

For a complete variant calling workflow with preprocessing, alignment, variant calling, and annotation:

```bash
nextflow run iarcbioinfo/pavoa-nf -profile apptainer \
  --input_file input.txt \
  --ref hg38.fasta \
  --trim \
  --bqsr \
  --known_sites dbsnp_138.hg38.vcf.gz \
  --known_sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --annovarDBlist hg38_listAVDB.txt
```

For dual index sequencing (UDseq)

```bash
nextflow run iarcbioinfo/pavoa-nf -profile apptainer \
  --input_file input.txt \
  --ref hg38.fasta \
  --trim \
  --umi \
  --known_sites dbsnp_138.hg38.vcf.gz \
  --known_sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --annovarDBlist hg38_listAVDB.txt
```

### Run only calling

The pair.txt file is a tabular file with three columns : SD (sample ID), normal (BAM file path), tumor (BAM file path).

```bash
nextflow run IARCbioinfo/pavoa-nf -entry dupcaller -profile apptainer \
  --input_file pairs.txt \
  --ref hg38.fasta \
  --annovarDBlist hg38_listAVDB.txt
```

## Output

| Type | Description |
|------|-------------|
| /    | Annotated variants in vcf and tabular format |
| BAM/ | folder with BAM and BAI files of alignments |
| VCF/ | folder with VCF files containing raw called variants |
| QC/BAM/ | alignment quality control reports |
| QC/duplicates/ | variant calling quality control reports |
| QC/multiqc_report.html | comprehensive MultiQC report |
| nf-pipeline_info/ | log files from all pipeline steps |


```
output/
├── BAM/
├── VCF/
├── QC/
│   ├── BAM/
│   ├── duplicates/
│   ├── multiqc_report.html
│   └── multiqc_data/
└── nf-pipeline_info/
```

## Workflow

The pavoa-nf pipeline performs the following major steps:

1. **Preprocessing (P)**
   - UMI trimming (optional, dupcaller trim)
   - Adapter trimming (optional, AdapterRemoval)
   - Quality control of reads (FastQC)
   
2. **Alignment (A)**
   - Read alignment to reference genome (BWA-MEM2)
   - Duplicate marking (samblaster)
   - Sorting and indexing (sambamba)
   - Base quality score recalibration (optional, GATK)
   - Alignment quality control (Qualimap)

3. **Variant calling (V)**
   - Variant calling using selected caller:
     - Dupcaller call (somatic variants, UDseq)

4. **Optimization (O)**

5. **Annotation (A)**
   - Variant annotation (Annovar)
   - Context annotation (Gama annot)

6. **Quality Control and Reporting**
   - Comprehensive MultiQC report

## Troubleshooting

### Common issues

### Resource requirements

- Minimum recommended resources: 8 CPUs, 32GB RAM
- For whole genome data: 16+ CPUs, 64GB+ RAM
- Temporary disk space: ~3x input file size

## Citation

If you use this pipeline, please cite:

- The pipeline: `pavoa-nf: Pipeline for variant detection. Preprocessing Alignment, Variant calling, Optimization, Annotation. https://github.com/IARCbioinfo/pavoa-nf`
- Nextflow: `Paolo Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017). doi:10.1038/nbt.3820`

Please also cite the individual tools used by the pipeline.

## Contributions

| Name | Email | Description |
|------|-------|-------------|
| Cahais Vincent* | cahaisv@iarc.who.int | Developer to contact for support |

## References

