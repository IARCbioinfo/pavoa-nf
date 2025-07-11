#! /usr/bin/env Rscript
#####################################################################################
#
# Title		: annovar_annot.r
# Author	: CahaisV@iarc.fr
# Date		: 03/07/2019
# Last Update	: 29cd /02/2024
#
#####################################################################################

library(data.table)
library(GetoptLong)

out=""
annovarDBpath="/data/databases/annovar/hg38db"
annovarBinPath="~/bin/annovar/"
threads=1
PASS="PASS"
tags=".bam|.vcf|.tsv|.somatic|.germline|.calls|.alt|.gz|.PASS|.filtered|.snvs|.indels|.fixed"

GetoptLong(matrix(c("input|i=s",          "vcf input file",	
                    "annovarDBlist|l=s",  "txt file listing annovar databases for annotation",
                    "annovarDBpath|a=s",  "path to annovarDB",
                    "annovarBinPath|b=s", "path to table_annovar.pl",
                    "threads|t=i",        "threads",
                    "PASS|p=s",           "filter tag",
                    "keepNonStandard|k",  "keep Non Standard Chromosomes",
                    "tags|g=s",           "strings paterns to remove in input file name"
), ncol=2, byrow=TRUE))


################
# Check parameters

print("Check input files")
#check annovarDBpath
if(!file.exists(annovarDBpath)){ stop( paste0("ERROR : ", annovarDBpath, " do not exit !") ) }
annovarDB=basename(annovarDBpath)
reference=gsub("db","",annovarDB)
print(paste0("Reference : ",reference ))

#check annovarDB list and create protocols and operation list
if(!file.exists(annovarDBlist)){ stop( paste0("ERROR : ", annovarDBlist, " do not exit !") ) }
dblist<-read.table(annovarDBlist,comment.char="#")
dblist$do<-apply(dblist, 1, function(row) file.exists(paste0(annovarDBpath,"/",row[1])) )
for( i in 1:nrow(dblist) ){ ifelse(dblist[i,3],print(paste0("VALID : ",dblist[i,1])),print(paste0("NOT FOUND : ", dblist[i,1])) ) }
dblist<-dblist[which(dblist$do==TRUE),]
protocols<-sub( "\\.txt", "",  sub(paste0(reference,"_"), "", dblist$V1) )
protocols<-paste(protocols,collapse = ",")
operations<-paste(dblist$V2,collapse = ",")
if(protocols==""){ stop( paste0("ERROR : there is no valid file in ", annovarDBlist, " !") ) }

# Check if input file exist and read it
if(!file.exists(input)){ stop( paste0("ERROR : ", input, " do not exit !") ) }

################
#Functions

filterVCF<-function(vcffile,PASS,keepNonStandard){
  
  PASS<-unlist(strsplit(PASS,","))

  # read vcf file
  vcf_lines <- readLines(vcffile)
  header_lines <- vcf_lines[startsWith(vcf_lines, "##")]
  data_lines   <- fread(vcffile,skip="#CHROM")[FILTER %in% PASS,]
  if(nrow(data_lines)==0){ stop("ERROR : None of the variants pass filters") }

  #keep only standard chromosomes
  if (!keepNonStandard){ 
    data_lines<-data_lines[ !grep("_",`#CHROM`),]
    data_lines<-data_lines[ !grep("M",`#CHROM`),]
  }

  #report number of chromosomes
  chrom=unique(data_lines$`#CHROM`)
  print(paste( "Chromosomes in vcf file = ", paste0(chrom,collapse=" ") ))

  writeLines(header_lines, "tmp.vcf")
  write.table(data_lines, "tmp.vcf", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  #compress
  err<-system("gzip -f tmp.vcf", ignore.stderr=T)
  if(err>0){stop(err)}

  return( colnames(data_lines) )
}

# restaure colnames in table_annovar ouputs
restaureColnames<-function(colNames,dir="./"){

  # multianno.vcf
  multianno <- list.files(path = dir, pattern="multianno.vcf")
  vcf_lines <- readLines(multianno)
  header_lines <- vcf_lines[startsWith(vcf_lines, "##")]
  data_lines   <- vcf_lines[!startsWith(vcf_lines, "##")]
  writeLines(c(header_lines,paste0(colNames,collapse="\t"),data_lines), multianno)

  # multianno.txt
  multianno <- list.files(path = dir, pattern="multianno.txt")
  txt_lines <-fread(multianno)
  col <- gsub("#","",colNames)
  colnames(txt_lines)[(ncol(txt_lines)-length(colNames)+1):ncol(txt_lines)] <- col
  fwrite(txt_lines, file = multianno, sep="\t")

  # avinput
  annotNames <- c("Chr","Start","End","Ref","Alt","s1","s2","s3")
  avinput  <-list.files(path = dir, pattern="avinput")
  data_lines <- readLines(avinput)
  col <- paste0(c(annotNames,col),collapse="\t")
  writeLines(c(col,data_lines), avinput)

}

###################
# MAIN

# Filter vcf
print("Filter VCF ...")
colNames<-filterVCF(input,PASS,keepNonStandard)

# Run annovar
print("Running Annovar ...")
annovarBin<-paste0( gsub("table_annovar.pl","",annovarBinPath), "/table_annovar.pl")
params=paste(" --buildver", reference , "--thread", threads, "--vcfinput --onetranscript --remove --otherinfo --protocol", protocols, "-operation", operations, "tmp.vcf.gz", annovarDBpath , sep=" ")
print("Run table_annovar.pl")
print(paste0(annovarBin,params))
err<-system2(command=annovarBin,args=params, stdout="annovar.stdout", stderr="annovar.stderr")

# PostProcess
print("Restaure columns names ...")
restaureColnames(colNames)

#rename files
tag <- gsub(tags,"",input)
outfiles<-list.files("./",pattern="tmp")
file.rename(outfiles,gsub("tmp",tag,outfiles))

# The end
print("The end !")
