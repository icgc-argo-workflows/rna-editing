# Credits -----------------------------------------------------------------
# Copyright (c) Center for Cancer Bioinformatics, Peking University Cancer Hosptial, All Rights Reserved.
# Author: Yang Yang, yangyang@bjcancer.org
# Create date: 2 Dec 2018.
# Last update date: 13 Dec 2018.

# HEADER ------------------------------------------------------------------
rm(list=ls())
options(stringsAsFactors=FALSE)
options(java.parameters ="-Xmx32384m")
ARG_MODE <- TRUE

if(ARG_MODE){
  args<-commandArgs(T)
  if (length(args)==0){
    cat("Usage: Rscript samtools_anno_summary.R [in_vcf_file] [analysisid] [out_file]","\n")
    cat("\t\tin_vcf_file: raw variant file from samtools","\n")
    cat("\t\tanalysisid: analysis id of sample","\n")
    cat("\t\tout_file: ","\n")
    q()
  }  
  in_vcf_file = args[1]
  analysisid = args[2]
  out_path = args[3]
  REDIportal_file = args[4]
  rmsk_region_file = args[5]
}else{
  in_vcf_file = "/mnt/dellfs/projects/RNA_editing/analysis_inhouse/samtools/W070_T2P37_3/W070_T2P37_3_filtered.vcf"
  analysisid = "W070_T2P37_3"
  out_path = "/mnt/dellfs/projects/RNA_editing/analysis_inhouse/samtools/W070_T2P37_3"
  REDIportal_file = "/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/rediportal/known_rnaediting_site.txt"
  rmsk_region_file = "/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/rmsk/rmsk_region_v2.gtf"
  #ref_gtf_file = ""
}

# Functions -------------------------------------------------------
source("/mnt/dellfs/pub/tools/ccb/ccb.helper.R")

# Input/Output Settings ---------------------------------------------------
# setwd(out_path)
## vcf file 
## in_vcf = read.file(in_vcf_file,header = F)
tmp_vcf <- readLines(in_vcf_file)
in_vcf <- read.file(in_vcf_file,header = F, stringsAsFactors = FALSE)
## filter for the columns names
tmp_vcf <- tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names <- unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
vcf_names[10] = analysisid
names(in_vcf) <- vcf_names

anno_vcf = in_vcf
anno_vcf$FILTER = ""
index = paste0(anno_vcf$`#CHROM`,"_",anno_vcf$POS)

# Annotation --------------------------------------------------------------
## Alu element annotation
rmsk = read.file(rmsk_region_file,header=F)
names(rmsk) = c("chromosome","start","end","strand","type","type_detail")
rmsk_alu = rmsk[grep("Alu",rmsk$type_detail),]
## select Alu
#rmsk=rmsk[rmsk$type %in% c("Simple_repeat","Satellite","Low_complexity"),]
for(i in c(1:nrow(anno_vcf))){
  if(nrow(rmsk_alu[rmsk_alu$chromosome == anno_vcf$`#CHROM`[i] & rmsk_alu$start<anno_vcf$POS[i] & rmsk_alu$end>anno_vcf$POS[i],])>0){
    anno_vcf$FILTER[i] = paste0(anno_vcf$FILTER[i],"ALU;")
  }
}
## Known RNA-editing
REDIportal = read.file(REDIportal_file,header=F)
names(REDIportal) = c("chromosome","position")
known_editing_pos = paste0(REDIportal$chromosome,"_",REDIportal$position)
if(nrow(anno_vcf[index %in% known_editing_pos,])>0){
  anno_vcf[index %in% known_editing_pos,]$FILTER=paste0(anno_vcf[index %in% known_editing_pos,]$FILTER,"KNOWN_EDITING;")
}
## A to I annotation
anno_vcf[anno_vcf$REF=="A" & anno_vcf$ALT=="G",]$FILTER=paste0(anno_vcf[anno_vcf$REF=="A" & anno_vcf$ALT=="G",]$FILTER,"A_to_I;")
anno_vcf[anno_vcf$REF=="T" & anno_vcf$ALT=="C",]$FILTER=paste0(anno_vcf[anno_vcf$REF=="T" & anno_vcf$ALT=="C",]$FILTER,"A_to_I;")
anno_vcf[anno_vcf$REF=="C" & anno_vcf$ALT=="T",]$FILTER=paste0(anno_vcf[anno_vcf$REF=="C" & anno_vcf$ALT=="T",]$FILTER,"C_to_U;")
anno_vcf[anno_vcf$REF=="G" & anno_vcf$ALT=="A",]$FILTER=paste0(anno_vcf[anno_vcf$REF=="G" & anno_vcf$ALT=="A",]$FILTER,"C_to_U;")
## output annotation vcf file
write.table(anno_vcf, file=paste0(out_path,"/",analysisid,"_filtered_anno.vcf"), quote=F, col.names = TRUE, row.names = FALSE, sep="\t")

# Summary -----------------------------------------------------------------
total=nrow(anno_vcf)
A_to_I=length(grep("A_to_I",anno_vcf$FILTER))
C_to_U=length(grep("C_to_U",anno_vcf$FILTER))
total_Alu=length(grep("ALU",anno_vcf$FILTER))
total_REDItools=length(grep("KNOWN_EDITING",anno_vcf$FILTER))
A_to_I_Alu=length(grep("ALU.*A_to_I",anno_vcf$FILTER))
A_I_REDIportal=length(grep("KNOWN_EDITING.*A_to_I",anno_vcf$FILTER))

summary=cbind(analysisid,total, A_to_I, C_to_U, total_Alu, total_REDItools, A_to_I_Alu, A_I_REDIportal)
##_filtering_annotation_summary.txt
write.table(summary, file=paste0(out_path,"/",analysisid,"_filtered_anno_summary.txt"), quote=F, col.names = TRUE, row.names = FALSE, sep="\t")

