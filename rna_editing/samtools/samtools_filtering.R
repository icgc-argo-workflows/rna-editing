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
    cat("Usage: Rscript samtools_filtering.R [in_vcf_file] [analysisid] [out_path] [somatic_maf_file]","\n")
    cat("\t\tin_vcf_file: raw variant file from samtools","\n")
    cat("\t\tanalysisid: analysis id of sample","\n")
    cat("\t\tout_path: ","\n")
    cat("\t\tsomatic_maf_file: TCGA/other","\n")
    q()
  }  
  in_vcf_file = args[1]
  analysisid = args[2]
  out_path = args[3]
  somatic_maf_file = args[4]
}else{
  in_vcf_file = "/mnt/dellfs/projects/RNA_editing/analysis_inhouse/samtools/W066T01O0/W066T01O0.vcf"
  analysisid = "W066T01O0"
  out_path = "/mnt/dellfs/projects/RNA_editing/analysis_inhouse/samtools/W066T01O0"
  somatic_maf_file = "/mnt/dellfs/projects/gc_multiregion/DNA/analysis_updating/cohort/maf/cohort.final.maf"
  # /mnt/dellfs/projects/RNA_editing/analysis_inhouse/samtools/W070_T2P37_2/W070_T2P37_2.vcf
  # snp_datasets_file = "/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/RNA-editing_combined_dbsnp.txt"
  # splicesites_file = "/mnt/dellfs/projects/RNA_editing/analysis_test/script/REDItools/resource/Gencode_annotation/gencode.v30lift37.splicesites.txt"
}

# common files
snp_position_file="/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/snp_position.txt"
splicesites_pos_file="/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/gencode.v30lift37.splicesites.pos.txt"
rmsk_region_file = "/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/rmsk_region.gtf"

# Functions -------------------------------------------------------
source("/mnt/dellfs/pub/tools/ccb/ccb.helper.R")

# Input/Output Settings ---------------------------------------------------
setwd(out_path)

# prepare files
# snp_datasets = read.file(snp_datasets_file,header=F)
# names(snp_datasets) = c("chromosome","position","ref","alt")
# snp_datasets$chromosome=gsub("chr", "", snp_datasets$chromosome) 
# snp_position=paste0(snp_datasets$chromosome,"_",snp_datasets$position)
# snp_position_file="/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/snp_position.txt"
# write.table(snp_position,file=file, quote=F, col.names = FALSE, row.names = FALSE, sep="\t",append = F)
# splicesites = read.table(splicesites_file,header=F,sep=" ")
# names(splicesites) = c("chromosome","start","end","A_D","strand")
# splicesites_pos_file="/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/gencode.v30lift37.splicesites.pos.txt"
# for(i in c(1:nrow(splicesites))){
#   chr = splicesites$chromosome[i]
#   start = splicesites$start[i]-3
#   end = splicesites$end[i]+3
#   tmp = cbind(rep(chr,length(end-start+1)),c(start:end))
#   write.table(tmp, file=splicesites_pos_file, quote=F, col.names = FALSE, row.names = FALSE, sep="\t",append = TRUE)
# }
# rmsk_region_file = "/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/rmsk_region.gtf"

## somatic file
if(somatic_maf_file == "TCGA"){
  somatic_file="/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/mc3_somatic_snp_cufoff2.txt"
  somatic_snp=read.file(somatic_file,header=F)
  names(somatic_snp)=c("chromosome","start","end")
}else{
  somatic_maf = read.file(somatic_maf_file,header=T)
  somatic_snp = subset(somatic_maf,Variant_Type=="SNP" & Tumor_Sample_Barcode==analysisid,select=c("Chromosome","Start_Position","End_Position"))
  names(somatic_snp)=c("chromosome","start","end")
}

## vcf file 
tmp_vcf<-readLines(in_vcf_file)
in_vcf<-read.table(in_vcf_file, header=F,stringsAsFactors = FALSE)
### filter for the columns names
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
vcf_names[10]=analysisid
names(in_vcf)<-vcf_names

# the candidate variant sites should be with base-quality ≥ 20, mapping quality ≥ 50, mapped reads ≥ 4, variant-supporting reads ≥ 3, 
# and mismatch frequencies (variant-supporting-reads/mapped-reads) ≥ 0.1.
# remove INDEL and homozygous loci;
filtering_quality=function(in_vcf, basequality=20, mapping_quality=20, mapped_reads=4, variant_supporting_reads=3, mismatch_frequencies=0.1){
  ## remove INDEL and homozygous loci
  vcf=in_vcf[-grep("INDEL",in_vcf$INFO),]
  vcf=vcf[-grep("1/1",vcf[,10]),]
  ## filtering QUAL
  vcf$FILTER=""
  vcf[vcf$QUAL<basequality,]$FILTER = "QUAL;"
  for(i in c(1:nrow(vcf))){
    FILTER=""
    INFO=strsplit(vcf$INFO[i],split=";",fixed = TRUE)[[1]]
    ## filtering mapping quality ≥ 50
    MQ=grep("MQ=",INFO,value=T)
    if(as.numeric(strsplit(MQ,split="=",fixed = TRUE)[[1]][2]) < mapping_quality){
      FILTER=paste0(FILTER,"MQ;")
    }
    ## filtering mapped reads ≥ 4
    DP=grep("DP=",INFO,value=T)
    if(as.numeric(strsplit(DP,split="=",fixed = TRUE)[[1]][2]) < mapped_reads){
      FILTER=paste0(FILTER,"DP;")
    }
    ## variant-supporting reads ≥ 3
    DP4 = grep("DP4=",INFO,value=T)
    DP4 = strsplit(strsplit(DP4,split="=",fixed = TRUE)[[1]][2],split=",",fixed = TRUE)[[1]]
    ref_reads = as.numeric(DP4[1])+as.numeric(DP4[2])
    alt_reads = as.numeric(DP4[3])+as.numeric(DP4[4])
    if(alt_reads < variant_supporting_reads){
      FILTER = paste0(FILTER,"DPALT;")
    }
    ## mismatch frequencies (variant-supporting-reads/mapped-reads) ≥ 0.1
    if(alt_reads/(ref_reads+alt_reads)<mismatch_frequencies){
      FILTER = paste0(FILTER,"ALTFREQ;")
    }
    vcf$FILTER[i] = paste0(vcf$FILTER[i],FILTER)
    # add
    # Statistical tests based on the binomial distribution B(n, p). p denotes the background mismatch rate of each transcriptome sequencing, 
    # and n denotes sequencing depth on this site.     
  }
  return(vcf)
  #vcf[vcf$FILTER=="",]$FILTER=="PASS"
}

anno_vcf = filtering_quality(in_vcf)
write.table(anno_vcf, file = paste0(out_path,"/",analysisid,"_quality_anno.vcf"), quote=F, col.names = TRUE, row.names = FALSE, sep="\t")
anno_vcf = anno_vcf[anno_vcf$FILTER=="",]
quality_filtering=nrow(in_vcf)-nrow(anno_vcf)

# Discard the sites present in combined DNA SNP datasets (dbSNP v.138, 1000 Genome SNP phase 3, human Dutch populations, 
# and BGI in-house data; combined datasets deposited at: ftp://ftp.genomics.org.cn/pub/icgc-pcawg3).
snp_position = read.file(snp_position_file,header=F)
snp_position = snp_position$V1
index = paste0(anno_vcf$`#CHROM`,"_",anno_vcf$POS)
#index_t = index[index %in% snp_position]
anno_vcf[index %in% snp_position,]$FILTER = paste0(anno_vcf[index %in% snp_position,]$FILTER,"dbSNP;")

# Estimate strand bias and filter out variants with strand bias based on two-tailed Fisher’s exact test. 

# Estimate and filter out variants with position bias, such as sites only found at the 3′ end or at 5′ end of a read.

# Discard the variation site in simple repeat region or homopolymer region or <5 bp from splicing site. 
rmsk = read.file(rmsk_region_file,header=F)
names(rmsk) = c("chromosome","HGNC","type","start","end")
## exclude simple repeat/low_complexity/Satellite
rmsk=rmsk[rmsk$type %in% c("Simple_repeat","Satellite","Low_complexity"),]
for(i in c(1:nrow(anno_vcf))){
  if(nrow(rmsk[rmsk$chromosome == anno_vcf$`#CHROM`[i] & rmsk$start<anno_vcf$POS[i] & rmsk$end>anno_vcf$POS[i],])>0){
    anno_vcf$FILTER[i] = paste0(anno_vcf$FILTER[i],"RMSK;")
  }
}

splicesites_pos = read.file(splicesites_pos_file,header=F)
names(splicesites_pos) = c("chromosome","position")
splicesites_pos = unique(splicesites_pos)
splicesites_pos_index = paste0(splicesites_pos$chromosome,"_",splicesites_pos$position)
anno_vcf[index %in% splicesites_pos_index,]$FILTER = paste0(anno_vcf[index %in% splicesites_pos_index,]$FILTER,"SPLICE;")

# retain a candidate variant site if at least 90% of its variant-supporting reads are realigned to this site. 
# Finally, all high confident RNA-editing sites were annotated by ANNOVAR.

# To remove the possibility of an RNA-editing variant being a somatic variant, the variant sites are positionally 
# filtered against PCAWG WGS somatic variant calls
# cutoff: >=2 samples detected.
## set cutoff for somatic mutation
somatic_snp_pos=paste0(gsub("[[:space:]]", "", somatic_snp$chromosome),"_",gsub("[[:space:]]", "", somatic_snp$start))
#somatic_snp_pos=paste0(gsub("[[:space:]]", "", somatic_snp$chromosome),"_",somatic_snp$start-1)
if(nrow(anno_vcf[index %in% somatic_snp_pos,])>0){
  anno_vcf[index %in% somatic_snp_pos,]$FILTER=paste0(anno_vcf[index %in% somatic_snp_pos,]$FILTER,"SOMATIC;")
}

write.table(anno_vcf, file = paste0(out_path,"/",analysisid,"_position_anno.vcf"), quote=F, col.names = TRUE, row.names = FALSE, sep="\t")

# Final vcf ---------------------------------------------------------------
filtered_vcf=anno_vcf[anno_vcf$FILTER=="",]
write.table(filtered_vcf, file=paste0(out_path,"/",analysisid,"_filtered.vcf"), quote=F, col.names = TRUE, row.names = FALSE, sep="\t")

# Summary filtering -------------------------------------------------------
raw=nrow(in_vcf)
RMSK_filtering=length(grep("RMSK",anno_vcf$FILTER))
SPLICE_filtering=length(grep("SPLICE",anno_vcf$FILTER))
dbSNP_filtering=length(grep("dbSNP",anno_vcf$FILTER))
SOMATIC_filtering=length(grep("SOMATIC",anno_vcf$FILTER))
final=length(filtered_vcf$POS)
summary=cbind(analysisid,raw, quality_filtering, RMSK_filtering, SPLICE_filtering, dbSNP_filtering, SOMATIC_filtering, final)
write.table(summary, file=paste0(out_path,"/",analysisid,"_filtering_summary.txt"), quote=F, col.names = TRUE, row.names = FALSE, sep="\t")
