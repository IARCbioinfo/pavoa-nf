#!/usr/bin/env Rscript

# Variant Filtering Pipeline - Batch Processing
# Usage: Rscript filter_variants.R [analyse_folder] [cov_n] [cov_t] [min_vaf_t] [max_vaf_t] [min_cov_alt_t] [max_cov_alt_t]

library(readr)
library(dplyr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parameters with defaults
sample <- if (length(args) >= 1) args[1] else "analyse"
cov_n_thresh <- if (length(args) >= 2) as.numeric(args[2]) else 10
cov_t_thresh <- if (length(args) >= 3) as.numeric(args[3]) else 10
min_vaf_t_thresh <- if (length(args) >= 4) as.numeric(args[4]) else 0.1
max_vaf_t_thresh <- if (length(args) >= 5) as.numeric(args[5]) else 1
min_cov_alt_t_thresh <- if (length(args) >= 6) as.numeric(args[6]) else 3
max_cov_alt_t_thresh <- if (length(args) >= 7) as.numeric(args[7]) else 1000

cat("Processing sample:", sample, "\n")
cat("Thresholds - Cov_N:", cov_n_thresh, "Cov_T:", cov_t_thresh,
    "min VAF_T:", min_vaf_t_thresh, "max VAF_T:", max_vaf_t_thresh,
    "min Cov_alt_T:", min_cov_alt_t_thresh,
    "max Cov_alt_T:", max_cov_alt_t_thresh, "\n")

# Create output directory
filtered_dir <- file.path("filtered.1")
if (!dir.exists(filtered_dir)) {
  dir.create(filtered_dir, recursive = TRUE)
}

# Write filtering parameters to log file
log_file <- file.path(filtered_dir, "filtering.log")
log_lines <- c(
  paste0("Cov_N >= ", cov_n_thresh),
  paste0("Cov_T >= ", cov_t_thresh),
  paste0("min VAF_T >= ", min_vaf_t_thresh),
  paste0("max VAF_T <= ", max_vaf_t_thresh),
  paste0("min_Cov_alt_T >= ", min_cov_alt_t_thresh),
  paste0("max_Cov_alt_T <= ", max_cov_alt_t_thresh),
  paste0("Cov_alt_N <= 0"),
  paste0("FILTER: PASS"),
  paste0("is_SNP == 0"),
  "blacklist.v2: . (if present)",
  "rmsk: . (if present)",
  "simpleRepeat: . (if present)"
)
writeLines(log_lines, log_file)


# Function to filter and create bed file
filter_variants <- function(tsv_file, vcf_file) {
  cat("\nProcessing:", basename(tsv_file), "\n")

  vcf_out <- file.path(paste0(filtered_dir, "/", basename(vcf_file) ))
  vcf_out <- gsub("\\.gz$", "", vcf_out)  # Remove .gz if present for output naming

  # Step 1: Filter TSV
  filtered <- read_tsv(tsv_file, show_col_types = FALSE) %>%
    filter(Cov_N >= cov_n_thresh,
           Cov_T >= cov_t_thresh,
           VAF_T >= min_vaf_t_thresh,
           VAF_T <  max_vaf_t_thresh,
           Cov_alt_T >= min_cov_alt_t_thresh,
           Cov_alt_T <= max_cov_alt_t_thresh,
           Cov_alt_N <= 0,
           is_SNP == 0,
           FILTER == "PASS"
    )

  if ("blacklist.v2" %in% names(filtered)) {
    filtered <- filtered %>% filter(blacklist.v2 == ".")
  }

  if ("rmsk" %in% names(filtered)) {
    filtered <- filtered %>% filter(rmsk == ".")
  }

  if ("simpleRepeat" %in% names(filtered)) {
    filtered <- filtered %>% filter(simpleRepeat == ".")
  }

  cat("Filtered variants:", nrow(filtered), "\n")
  
  if (nrow(filtered) == 0) {
    cat("No variants passed filters\n")
    file.create(vcf_out)  # Create empty VCF file
    return()
  }

  # Step 2: Create BED file
  bed_file <- file.path(gsub(vcf_file, pattern = "\\.vcf$|*.vcf.gz$", replacement = ".bed"))
  bed_data <- filtered %>% select(1, 2, 3)  # First 3 columns: chr, start, end
  write_tsv(bed_data, bed_file, col_names = FALSE)

  # Step 3: Filter VCF using bedtools
  cmd <- sprintf("bedtools intersect -a %s -b %s -header > %s", vcf_file, bed_file, vcf_out)
  system(cmd)

  cat("Created:", basename(vcf_out), "\n")
}

files<-list.files("./", pattern=sample, full.names = TRUE)
vcf <- files[grep("vcf$|vcf.gz$", files)]
vcf <- grep("multianno", vcf, invert = TRUE, value = TRUE)[1]
tsv <- files[grep("tsv$", files)][1]
#tsv <- list.files("./", pattern = "*.1.tsv$", full.names = TRUE)[1]
#vcf <- list.files("./", pattern = "*.vcf$|*.vcf.gz$", full.names = TRUE)
#vcf <- grep("multianno", vcf, invert = TRUE, value = TRUE)[1]

filter_variants(tsv[1], vcf[1])

cat("\n=== Pipeline completed ===\n")
cat("Results saved in:", filtered_dir, "\n")
