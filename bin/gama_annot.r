#! /usr/bin/env Rscript
#####################################################################################
#
# Title		: gama_annot.r
# Author	: CahaisV@iarc.fr
# Date		: 31/07/2020
# Last Update	: 29/02/2024
#
#####################################################################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GetoptLong))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))

annovarDBpath <- "/data/databases/annovar/mm10db"

GetoptLong(matrix(c("annovarDBpath|a=s", "path to annovarDB"), ncol = 2, byrow = TRUE))

annovarDB <- basename(annovarDBpath)
reference <- gsub("db", "", annovarDB)

##########################################
# order_variants
order_variants<-function(tab){

  print("order_variants")
  # order lines
  tab %>% mutate(num = Chr) %>%
    mutate(num = ifelse(num == "chrX", "1023", num)) %>%
    mutate(num = ifelse(num == "chrY", "1024", num)) %>%
    mutate(num = ifelse(num == "chrM", "1025", num)) %>%
    mutate(num = as.integer(gsub("chr", "", num))) %>%
    arrange(num, Start, End, Ref, Alt) %>%
    filter(!is.na(Chr), Chr != "NA") %>%
    select(-num)

}

##########################################
# get_caller_name
get_caller_name <- function(headers) {

  print("get_caller_name")
  callers=c("strelka", "Mutect2", "octopus", "needlestack","dupcaller")
  for (caller in callers) {
    if (any(grepl(caller, headers))) {
      return(caller)
    }
  }
  return("dupcaller") # dupcaller do not sign its vcf :-(
}

##########################################
# set_caller_name
set_caller_name<-function(tab,dir="./"){

  print("set_caller_name")
  vcffile<-list.files(path = dir, pattern="multianno.vcf")
  vcffile <- readLines(vcffile)
  header_lines <- vcffile[startsWith(vcffile, "##")]
  tab %>% mutate(caller= get_caller_name(header_lines) )
}

##########################################
# Retrieve Strand from refGene file
# It cannot work if refGene is not used !
getStrand <- function(avtmp) {

  print("getStrand")
  if (is.null(avtmp$Gene.refGene)) {
    return(avtmp)
  }

  #Get strand information for each genes in refGene
  refGene <- fread(paste0(annovarDBpath, "/", reference, "_refGene.txt"))
  tmp <- unique(refGene[, c(13, 4)])
  colnames(tmp) <- c("symbol", "Strand")
  dup <- tmp$symbol[duplicated(tmp$symbol)]
  tmp[symbol %in% dup, ]$Strand <- "+/-"
  tmp <- tmp[!duplicated(tmp$symbol), ]

  # Add refGene strand info to avtmp
  setDT(avtmp)
  avtmp$symbol <- gsub("\\(.*", "", avtmp$Gene.refGene)
  avtmp[, symbol2 := tstrsplit(symbol, ";", keep = 1)]
  avtmp <- merge(avtmp, tmp, by = "symbol", all.x = T)
  avtmp[, symbol := NULL]
  avtmp[, symbol2 := NULL]
  return(as_tibble(avtmp))
}


##########################################
# Retrieve Context from fasta reference file
getContextAnnotation <- function(avtmp) {

  print("getContextAnnotation")
  reffile <- list.files(path = annovarDBpath, pattern = paste0(reference, ".fa"), full.names = T)
  ref <- readDNAStringSet(reffile)
  setDT(avtmp)
  avtmp[, context := getContext(ref, CHROM, POS, 10)]
  avtmp$context <- as.character(avtmp$context)
  avtmp[, trinucleotide_context := substr(context, 10, 12)]
  avtmp$trinucleotide_context <- sub("(.).(.)", "\\1x\\2", avtmp$trinucleotide_context)
  return(as_tibble(avtmp))
}

getContext <- function(ref, chr, pos, win) {
  return(mclapply(1:length(chr), function(x) toString(subseq(ref[[chr[x]]], max(1, as.numeric(pos[x]) - win), min(length(ref[[chr[x]]]), as.numeric(pos[x]) + win)))))
}

#######################################################
# Get VAF                                             #
#######################################################

# getVAF_Octopus
getVAF_Octopus <- function(vcf) {

  print("getVAF_Octopus")
  # Check format
  Format <- (vcf %>% extract(FORMAT, c("ADAFADP"), "([^:]+:[^:]+:[^:]+):FT") %>% pull(ADAFADP))[1]
  stopifnot(assertthat::are_equal(Format, "AD:AF:ADP"))

  vcf <- vcf %>%
    rename_with(.cols = ncol(.), ~"TUMOR") %>%
    rename_with(.cols = ncol(.) - 1, ~"NORMAL") %>%
    extract(TUMOR, c("AD", "AF", "Cov_T"), ":([^:]+):([^:]+):([^:]+):[^:]+$", remove = F) %>%
    separate(AD, into = c("Cov_ref_T", "Cov_alt_T"), extra = "drop", sep = ",") %>%
    dplyr::select(-Cov_ref_T) %>%
    separate(AF, into = c("VAF_ref_T", "VAF_T"), extra = "drop", sep = ",") %>%
    dplyr::select(-VAF_ref_T) %>%
    mutate(Cov_T = as.numeric(Cov_T), Cov_alt_T = as.numeric(Cov_alt_T), VAF_T = as.numeric(VAF_T), QUAL = as.numeric(QUAL)) %>%
    extract(NORMAL, c("AD", "AF", "Cov_N"), ":([^:]+):([^:]+):([^:]+):[^:]+$", remove = F) %>%
    separate(AD, into = c("Cov_ref_N", "Cov_alt_N"), extra = "drop", sep = ",") %>%
    dplyr::select(-Cov_ref_N) %>%
    separate(AF, into = c("VAF_ref_N", "VAF_N"), extra = "drop", sep = ",") %>%
    dplyr::select(-VAF_ref_N) %>%
    mutate(Cov_N = as.numeric(Cov_N), Cov_alt_N = as.numeric(Cov_alt_N), VAF_N = as.numeric(VAF_N), QUAL = as.numeric(QUAL)) %>%
    relocate("Cov_N", "Cov_T", "VAF_N", "VAF_T", "Cov_alt_T", "Cov_alt_N", .after = last_col())

  return(vcf)
}

# getVAF_Strelka
getVAF_Strelka <- function(vcf) {
  
  print("getVAF_Strelka")
  # Check format
  Format <- (vcf %>% pull(FORMAT))[1]
  if (assertthat::are_equal(Format, "GT:DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50")) {
    return(getVAF_StrelkaIndels(vcf))
  } else if (assertthat::are_equal(Format, "GT:DP:FDP:SDP:SUBDP:AU:CU:GU:TU")) {
    return(getVAF_StrelkaSNV(vcf))
  } else if (assertthat::are_equal(Format, "GT:DP:FGT:DP:SGT:DP:SUBGT:DP:AU:CU:GU:TU")) {
    return(getVAF_StrelkaSNV(vcf))
  } else {
    stop
  }
}

# getVAF_StrelkaSNV
getVAF_StrelkaSNV <- function(vcf) {

  print("getVAF_StrelkaSNV")
  # Check format
  # Format=( vcf %>% pull(FORMAT) )[1]
  # stopifnot( assertthat::are_equal(Format,"GT:DP:FDP:SDP:SUBDP:AU:CU:GU:TU") )

  vcf <- vcf %>%
    extract(TUMOR, c("Cov_T", "A", "C", "G", "T"), "[^:]+:([^:]+):.+:([^:]+):([^:]+):([^:]+):([^:])+$", remove = F) %>%
    mutate(Cov_alt_T = case_when(Alt == "A" ~ A, Alt == "C" ~ C, Alt == "G" ~ G, Alt == "T" ~ T, TRUE ~ ".")) %>%
    dplyr::select(-c(A, C, G, T)) %>%
    separate(Cov_alt_T, into = c("Cov_alt_T"), extra = "drop") %>%
    mutate(Cov_T = as.numeric(Cov_T), Cov_alt_T = as.numeric(Cov_alt_T), VAF_T = Cov_alt_T / Cov_T, QUAL = as.numeric(QUAL)) %>%
    extract(NORMAL, c("Cov_N", "A", "C", "G", "T"), "[^:]+:([^:]+):.+:([^:]+):([^:]+):([^:]+):([^:])+$", remove = F) %>%
    mutate(Cov_alt_N = case_when(Alt == "A" ~ A, Alt == "C" ~ C, Alt == "G" ~ G, Alt == "T" ~ T, TRUE ~ ".")) %>%
    dplyr::select(-c(A, C, G, T)) %>%
    separate(Cov_alt_N, into = c("Cov_alt_N"), extra = "drop") %>%
    mutate(Cov_N = as.numeric(Cov_N), Cov_alt_N = as.numeric(Cov_alt_N), VAF_N = Cov_alt_N / Cov_N, QUAL = as.numeric(QUAL)) %>%
    relocate("Cov_N", "Cov_T", "VAF_N", "VAF_T", "Cov_alt_T", "Cov_alt_N", .after = last_col())

  return(vcf)
}

# getVAF_StrelkaIndels
# VAF = TIR(1 tier) / DP
getVAF_StrelkaIndels <- function(vcf) {

  print("getVAF_StrelkaIndels")
  # Check format
  # Format=( vcf %>% pull(FORMAT) )[1]
  # stopifnot( assertthat::are_equal(Format,"GT:DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50") )

  vcf <- vcf %>%
    extract(TUMOR, c("Cov_T", "Cov_alt_T"), "[^:]+:([^:]+):[^:]+:[^:]+:([^:]+):", remove = F) %>%
    separate(Cov_alt_T, into = c("Cov_alt_T"), extra = "drop") %>%
    mutate(Cov_T = as.numeric(Cov_T), Cov_alt_T = as.numeric(Cov_alt_T), VAF_T = Cov_alt_T / Cov_T, QUAL = as.numeric(QUAL)) %>%
    extract(NORMAL, c("Cov_N", "Cov_alt_N"), "[^:]+:([^:]+):[^:]+:[^:]+:([^:]+):", remove = F) %>%
    separate(Cov_alt_N, into = c("Cov_alt_N"), extra = "drop") %>%
    mutate(Cov_N = as.numeric(Cov_N), Cov_alt_N = as.numeric(Cov_alt_N), VAF_N = Cov_alt_N / Cov_N, QUAL = as.numeric(QUAL)) %>%
    relocate("Cov_N", "Cov_T", "VAF_N", "VAF_T", "Cov_alt_T", "Cov_alt_N", .after = last_col())

  return(vcf)
}

# getVAF_Mutect
getVAF_Mutect <- function(vcf) {

  print("getVAF_Mutect")
  # Check format
  Format <- (vcf %>% extract(FORMAT, c("ADAF"), "[^:]+:([^:]+:[^:]+):") %>% pull(ADAF))[1]
  stopifnot(assertthat::are_equal(Format, "AD:AF"))

  cf <- which(colnames(vcf) == "FORMAT")
  colnames(vcf)[(cf + 1):(cf + 2)] <- c("TUMOR", "NORMAL") # this order is checked later

  vcf <- vcf %>%
    extract(TUMOR, c("Cov", "VAF_T"), "[^:]+:([^:]+):([^:]+):", remove = F) %>%
    separate(Cov, into = c("Cov_T", "Cov_alt_T"), extra = "drop") %>%
    mutate(Cov_T = as.numeric(Cov_T), Cov_alt_T = as.numeric(Cov_alt_T), VAF_T = as.numeric(VAF_T), QUAL = as.numeric(QUAL)) %>%
    mutate(Cov_T = Cov_T + Cov_alt_T) %>%
    extract(NORMAL, c("Cov", "VAF_N"), "[^:]+:([^:]+):([^:]+):", remove = F) %>%
    separate(Cov, into = c("Cov_N", "Cov_alt_N"), extra = "drop") %>%
    mutate(Cov_N = as.numeric(Cov_N), Cov_alt_N = as.numeric(Cov_alt_N), VAF_N = as.numeric(VAF_N), QUAL = as.numeric(QUAL)) %>%
    mutate(Cov_N = Cov_N + Cov_alt_N) %>%
    relocate("Cov_N", "Cov_T", "VAF_N", "VAF_T", "Cov_alt_T", "Cov_alt_N", .after = last_col())

  # Invert NORMAL and TUMOR when they are not in the correct order ( Cov_alt_T is suposed to be superior in TUMOR ! )
  if (sum(vcf$Cov_alt_T, na.rm=TRUE) < sum(vcf$Cov_alt_N,na.rm=TRUE)) {
    vcf <- vcf %>%
      dplyr::rename(TUMOR = NORMAL, NORMAL = TUMOR, Cov_T = Cov_N, Cov_N = Cov_T, VAF_T = VAF_N, VAF_N = VAF_T, Cov_alt_T = Cov_alt_N, Cov_alt_N = Cov_alt_T) %>%
      relocate("Cov_N", "Cov_T", "VAF_N", "VAF_T", "Cov_alt_T", "Cov_alt_N", .after = last_col())
  }

  return(vcf)
}


# getVAF_Needlestack
getVAF_Needlestack <- function(vcf) {

  print("getVAF_Needlestack")
  # Check format
  Format <- (vcf %>% extract(FORMAT, c("DPAOAF"), "GT:[^:]+:([^:]+:[^:]+:[^:]+:[^:]+):SB") %>% pull(DPAOAF))[1]
  stopifnot(assertthat::are_equal(Format, "DP:RO:AO:AF"))

  Format <- which(colnames(vcf) == "FORMAT")
  nbSamp <- (ncol(vcf) - Format) / 3

  # create sampleID base on position of sample columns
  sampID <- rep(seq(1, nrow(vcf) * nbSamp))
  for (i in 1:nbSamp) {
    v <- (i - 1) * nrow(vcf) + 1
    w <- (i - 1) * nrow(vcf) + nrow(vcf)
    sampID <- c(sampID, rep(seq(v, w), 2))
  }

  vcf <- vcf %>%
    gather("samples", "values", (Format + 1):ncol(.)) %>%
    mutate(sampID = sampID) %>%
    group_by(across(-all_of(c("samples", "values")))) %>%
    summarize(samples = paste(samples, collapse = " "), values = paste(values, collapse = " "), .groups = "drop") %>%
    separate(samples, into = c("Baseline", "Sample", "Sample_"), sep = " ") %>%
    separate(values, into = c("NORMAL", "TUMOR", "TUMOR_"), sep = " ") %>%
    extract(TUMOR, c("Cov_T", "Cov_alt_T", "VAF_T", "STATUS"), ":[^:]+:([^:]+):[^:]+:([^:]+):([^:]+):.+:([^:]+)$", remove = F) %>%
    mutate(Cov_T = as.numeric(Cov_T), Cov_alt_T = as.numeric(Cov_alt_T), VAF_T = as.numeric(VAF_T), QUAL = as.numeric(QUAL)) %>%
    extract(NORMAL, c("Cov_N", "Cov_alt_N", "VAF_N", "STATUS"), ":[^:]+:([^:]+):[^:]+:([^:]+):([^:]+):.+:([^:]+)$", remove = F) %>%
    mutate(Cov_N = as.numeric(Cov_N), Cov_alt_N = as.numeric(Cov_alt_N), VAF_N = as.numeric(VAF_N), QUAL = as.numeric(QUAL)) %>%
    relocate("Cov_N", "Cov_T", "VAF_N", "VAF_T", "Cov_alt_T", "Cov_alt_N", .after = last_col())

  return(vcf)
}


#getVAF_Dupcaller
getVAF_Dupcaller <- function(vcf){
  print("getVaf_Dupcaller")

  Format <- (vcf %>% pull(FORMAT))[1]
  stopifnot(assertthat::are_equal(Format, "GT:AC:RC:DP"))

  vcf <- vcf %>%
    extract(TUMOR, c("Cov_alt_T", "Cov_T"), "[^:]+:([^:]+):[^:]+:([^:]+)", remove = F) %>%
    mutate(Cov_T = as.numeric(Cov_T), Cov_alt_T = as.numeric(Cov_alt_T), VAF_T = Cov_alt_T / Cov_T) %>%
    extract(NORMAL, c("Cov_alt_N", "Cov_N"), "[^:]+:([^:]+):[^:]+:([^:]+)", remove = F) %>%
    mutate(Cov_N = as.numeric(Cov_N), Cov_alt_N = as.numeric(Cov_alt_N), VAF_N = Cov_alt_N / Cov_N) %>%
    relocate("Cov_N", "Cov_T", "VAF_N", "VAF_T", "Cov_alt_T", "Cov_alt_N", .after = last_col()) 

  return(vcf)
}


# get_VAF
get_VAF <- function(vcf) {

  print("getVAF")
  caller <- pull(vcf, caller)[1]
  switch(caller,
    strelka = return(getVAF_Strelka(vcf)),
    octopus = return(getVAF_Octopus(vcf)),
    needlestack = return(getVAF_Needlestack(vcf)),
    Mutect2 = return(getVAF_Mutect(vcf)),
    dupcaller = return(getVAF_Dupcaller(vcf)),
    stop("Unknown caller")
  )
}


###################################################################
# file_annot
###################################################################

# reverseComplement
reverseComplement <- function(change_strand, ctx) {
  revcomp <- list(
    AxA = "TxT", AxC = "GxT", AxG = "CxT", AxT = "AxT",
    CxA = "TxG", CxC = "GxG", CxG = "CxG", CxT = "AxG",
    GxA = "TxC", GxC = "GxC", GxG = "CxC", GxT = "AxC",
    TxA = "TxA", TxC = "GxA", TxG = "CxA", TxT = "AxA"
  )
  ifelse(change_strand, return(revcomp[[ctx]][[1]]), return(ctx))
}

reverseComplement_V <- Vectorize(reverseComplement)

# annotate_SNP
annotate_SNP <- function(dff) {
  # annotate human snp
  if ("ALL.sites.2015_08" %in% names(dff)) {
    dff <- dff %>% mutate(is_SNP = ifelse(ALL.sites.2015_08 >= 0.001 |
      ExAC_nontcga_ALL >= 0.001 |
      gnomAD_genome_ALL >= 0.001 |
      esp6500siv2_all >= 0.00001, 1, 0))
  }
  # annotate mouse snp
  if ("snp142" %in% names(dff)) {
    dff <- dff %>% mutate(is_SNP = ifelse(snp142 == ".", 0, 1))
  }

  # annotate rat rn6 snp
  if ("snp146" %in% names(dff)) {
    dff <- dff %>%
      mutate(is_SNP = ifelse(snp146 == ".", 0, 1)) %>%
      mutate(is_rpt = 0) %>%
      mutate(is_blacklist = 0) %>%
      mutate(is_hpoly = 0, is_str = 0)
  }

  # annotate rat rn7 snp
  if ("snp" %in% names(dff)) {
    dff <- dff %>%
      mutate(is_SNP = ifelse(snp == ".", 0, 1)) %>%
      mutate(is_rpt = ifelse(trf == ".", 0, 1)) %>%
      mutate(is_blacklist = ifelse(rmk == ".", 0, 1)) %>%
      mutate(is_hpoly = 0, is_str = 0)
  }

  return(dff)
}

# file_annot
file_annot <- function(df) {
  print("file_annot")
  df_var <- df %>%
    filter(!is.na(Chr)) %>%
    # mut_type
    mutate(mut_type = case_when(
      Ref == "-" | nchar(Alt) > 1 ~ "INS",
      Alt == "-" | nchar(Ref) > 1 ~ "DEL",
      TRUE ~ "SNV"
    )) %>%
    # mutation_class
    mutate(mut_class = case_when(
      (Ref == "C" & Alt == "A") | (Ref == "G" & Alt == "T") ~ "C>A",
      (Ref == "C" & Alt == "G") | (Ref == "G" & Alt == "C") ~ "C>G",
      (Ref == "C" & Alt == "T") | (Ref == "G" & Alt == "A") ~ "C>T",
      (Ref == "T" & Alt == "A") | (Ref == "A" & Alt == "T") ~ "T>A",
      (Ref == "T" & Alt == "C") | (Ref == "A" & Alt == "G") ~ "T>C",
      (Ref == "T" & Alt == "G") | (Ref == "A" & Alt == "C") ~ "T>G",
      TRUE ~ "-"
    )) %>%
    # add collapse context
    mutate(clps_context = as.character(reverseComplement_V(Ref %in% c("A", "G"), trinucleotide_context))) %>%
    # add collapse strand
    mutate(clps_tx_status = case_when(
      Ref %in% c("A", "G") & Strand == "-" ~ "UnTranscribed",
      Ref %in% c("A", "G") & Strand == "+" ~ "Transcribed",
      Strand == "+" ~ "UnTranscribed",
      Strand == "-" ~ "Transcribed",
      TRUE ~ "NA"
    )) %>%
    mutate(clps_strand = substr(clps_tx_status, 1, 1)) %>%
    # Annotate SNP
    annotate_SNP()

  return(df_var)
}

##################################
# main                           #
##################################

print("Load avinput")
multianno<-list.files(path = "./", pattern="multianno.txt")[1]

# set output file name
out <- gsub(".txt", ".1.tsv", multianno)
print(paste0("Output file name : ", out))

multianno <- multianno %>% read_tsv()
multianno <- multianno %>% order_variants()
multianno <- multianno %>% set_caller_name()
multianno <- multianno %>% dplyr::filter(ALT!=".") 
multianno <- multianno %>% getStrand()
multianno <- multianno %>% getContextAnnotation()
multianno <- multianno %>% get_VAF()
multianno <- multianno %>% file_annot()

print("Write output")
write_tsv(multianno, file = out)
print("The end !")
