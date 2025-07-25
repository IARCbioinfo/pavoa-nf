#!/usr/bin/env Rscript

# Variant Filtering Pipeline - Batch Processing
# Usage: Rscript filter_variants.R [analyse_folder] [cov_n] [cov_t] [vaf_t] [cov_alt_t]

library(readr)
library(dplyr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parameters with defaults
sample <- if (length(args) >= 1) args[1] else "analyse"
cov_n_thresh <- if (length(args) >= 2) as.numeric(args[2]) else 10
cov_t_thresh <- if (length(args) >= 3) as.numeric(args[3]) else 10
vaf_t_thresh <- if (length(args) >= 4) as.numeric(args[4]) else 0.1
cov_alt_t_thresh <- if (length(args) >= 5) as.numeric(args[5]) else 3

cat("Processing sample:", sample, "\n")
cat("Thresholds - Cov_N:", cov_n_thresh, "Cov_T:", cov_t_thresh, 
    "VAF_T:", vaf_t_thresh, "Cov_alt_T:", cov_alt_t_thresh, "\n")

# Create output directory
filtered_dir <- file.path(sample)
if (!dir.exists(filtered_dir)) {
  dir.create(filtered_dir, recursive = TRUE)
}

# Function to filter and create bed file
filter_variants <- function(tsv_file, vcf_file, sample) {
  cat("\nProcessing:", basename(tsv_file), "\n")
  
  # Step 1: Filter TSV
  filtered <- read_tsv(tsv_file, show_col_types = FALSE) %>%
    filter(avsnp150 == ".",
           ALL.sites.2015_08 == ".",
           AMR.sites.2015_08 == ".",
           Cov_N >= cov_n_thresh,
           Cov_T >= cov_t_thresh,
           VAF_T >= vaf_t_thresh,
           Cov_alt_T >= cov_alt_t_thresh,
           Cov_alt_N == 0,
           is_SNP == 0)
  
  cat("Filtered variants:", nrow(filtered), "\n")
  
  if (nrow(filtered) == 0) {
    cat("No variants passed filters\n")
    return()
  }
  
  # Step 2: Create BED file
  bed_file <- file.path(gsub(vcf_file, pattern = "\\.vcf$", replacement = ".bed"))
  bed_data <- filtered %>% select(1, 2, 3)  # First 3 columns: chr, start, end
  write_tsv(bed_data, bed_file, col_names = FALSE)
  
  # Step 3: Filter VCF using bedtools
  vcf_out <- file.path(paste0(sample,"/", sample, ".vcf"))
  cmd <- sprintf("bedtools intersect -a %s -b %s -header > %s", vcf_file, bed_file, vcf_out)
  system(cmd)
  
  cat("Created:", basename(vcf_out), "\n")
}

tsv<-list.files("./", pattern = ".*_multianno.1.tsv$", full.names = TRUE)
vcf<-list.files("./", pattern = ".*_multianno.vcf$", full.names = TRUE)

filter_variants(tsv[1], vcf[1], sample)

cat("\n=== Pipeline completed ===\n")
cat("Results saved in:", filtered_dir, "\n")
