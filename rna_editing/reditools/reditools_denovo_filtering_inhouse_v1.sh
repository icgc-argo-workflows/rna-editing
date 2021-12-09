#!/bin/bash
##############################################################################################
# Copyright(c) bjcancer, all rights reserved
# @file        : reditools_annotation.sh
# @author      : yangyang
# @email       : yangyang@bjcancer.org
# @revision    : 2021-07-8 03:14:20
# @description :       
###############################################################################################
if [ $# != 8 ] ; then 
    echo -e "Usage: reditools_denovo_filtering_inhouse.sh [reditools_path] out_path analysis_id somatic_maf rmsk_file snp151_file rediportal_file splicesites_gtf_file"
    echo -e "\t\treditools_path: para1"
    echo -e "\t\tout_path: para2"
    echo -e "\t\tanalysis_id: para3"
    echo -e "\t\tsomatic_maf: para4"
    echo -e "\t\trmsk_file: para5"
    echo -e "\t\tsnp151_file: para6"
    echo -e "\t\trediportal_file: para7"
    echo -e "\t\tsplicesites_gtf_file: para8"
    echo -e "example: \n"
exit 1;
fi

reditools_path=$1
out_path=$2
analysis_id=$3
somatic_maf=$4
rmsk_file=$5
snp151_file=$6
rediportal_file=$7
splicesites_gtf_file=$8

mkdir $out_path/${analysis_id}
result=$out_path/${analysis_id}/outTable_${analysis_id}.out

### maf to gtf
if [ $somatic_maf == "TCGA" ]; then
    somatic_gtf=/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/somatic/somatic.sorted.gtf.gz
else
    somatic_gtf=$out_path/${analysis_id}/somatic.sorted.gtf.gz
    cut -f5,6,7,8,10,16 $somatic_maf|grep -v '^#'|grep -w SNP|grep -w $analysis_id|awk '{print "chr"$1"\tsomatic\tsomatic\t"$2"\t"$3"\t.\t"$4"\t.\tgene_id \"somatic\"; transcript_id \"somatic\";"}' >$out_path/${analysis_id}/somatic.gtf
    sort -k1,1 -k4,4n $out_path/${analysis_id}/somatic.gtf > $out_path/${analysis_id}/somatic.sorted.gtf
    bgzip $out_path/${analysis_id}/somatic.sorted.gtf
    tabix -p gff $out_path/${analysis_id}/somatic.sorted.gtf.gz
fi

# Enter the REDItools denovo output folder:
in_file=`ls $reditools_path/denovo_${analysis_id}*/outTable_*`
# Exclude invariant positions and did not pass fisher test and genotyping
awk 'FS="\t" {if ($8!="-" && $10<0.05 && $9>=0.1) print}' $in_file >$result

# quality filtering: [base-quality ≥ 20, mapping quality ≥ 50], mapped reads ≥ 4, variant-supporting reads ≥ 3, mismatch frequencies (variant-supporting-reads/mapped-reads) ≥ 0.1
# remove INDEL and homozygous loci;
python /REDItools/accessory/selectPositions.py -i $result -c 4 -v 3 -f 0.1 -r -o ${result}.quality.sel

# database annotation and filtering: rmsk; snp; splicesite; SOMATIC of private maf; REDIportal annotation
# rmsk(11-12); snp (13-14); splicesite (15-16); somatic (17-18); REDIportal (19-20); 
python /REDItools/accessory/AnnotateTable.py -a $rmsk_file -n rmsk -r chr -i ${result}.quality.sel -o ${result}.quality.sel.rmsk -u
python /REDItools/accessory/AnnotateTable.py -a $snp151_file -n snp151 -r chr -i ${result}.quality.sel.rmsk -o ${result}.quality.sel.rmsk.snp -u
python /REDItools/accessory/AnnotateTable.py -a $splicesites_gtf_file -n splicesite -r chr -i ${result}.quality.sel.rmsk.snp -o ${result}.quality.sel.rmsk.snp.splicesite -u
python /REDItools/accessory/AnnotateTable.py -a $somatic_gtf -n somatic -r chr -i ${result}.quality.sel.rmsk.snp.splicesite -o ${result}.quality.sel.rmsk.snp.splicesite.somatic -u
python /REDItools/accessory/AnnotateTable.py -a $rediportal_file -n known_editing_site -r chr -i ${result}.quality.sel.rmsk.snp.splicesite.somatic -o ${result}.quality.sel.rmsk.snp.splicesite.somatic.rediportal -u

# final filtering output
awk 'FS="\t" {if ($1!="M" && $11!="Simple_repeat" && $11!="Low_complexity" && $13=="-" && $15=="-" && $17=="-") print}' ${result}.quality.sel.rmsk.snp.splicesite.somatic.rediportal  > ${result}.editing
awk 'FS="\t" {if (substr($12,1,3)=="Alu") print}' ${result}.editing > ${result}.editing.alu
awk 'FS="\t" {if (substr($12,1,3)!="Alu") print}' ${result}.editing > ${result}.editing.nonalu

# summary
echo -e "raw\tquality\trmsk\tsnp\tsplicesite\tsomatic\tfinal\tA_to_I\tC_to_U\ttotal_alu\ttotal_rediportal\tA_I_Alu\tA_I_REDIportal" >$out_path/${analysis_id}/${analysis_id}_filtering_summary.txt
raw=`wc -l $in_file|cut -d' ' -f1`
quality=`wc -l $result|cut -d' ' -f1`
rmsk=`cut -f 11 ${result}.quality.sel.rmsk.snp.splicesite.somatic.rediportal|grep 'Simple_repeat\|Low_complexity'|wc -l|cut -d' ' -f1`
snp=`cut -f 13 ${result}.quality.sel.rmsk.snp.splicesite.somatic.rediportal|grep 'snp'|wc -l|cut -d' ' -f1`
splicesite=`cut -f 15 ${result}.quality.sel.rmsk.snp.splicesite.somatic.rediportal|grep 'splicesite'|wc -l|cut -d' ' -f1`
somatic=`cut -f 15 ${result}.quality.sel.rmsk.snp.splicesite.somatic.rediportal|grep 'somatic'|wc -l|cut -d' ' -f1`
final=`wc -l ${result}.editing|cut -d' ' -f1`
A_to_I=`grep 'AG\|TC' ${result}.editing|wc -l|cut -d' ' -f1`
C_to_U=`grep 'GA\|CT' ${result}.editing|wc -l|cut -d' ' -f1`
total_alu=`wc -l ${result}.editing.alu|cut -d' ' -f1`
total_rediportal=`cut -f 19 ${result}.editing|grep 'ed'|wc -l|cut -d' ' -f1`
A_I_Alu=`grep 'AG\|TC' ${result}.editing.alu|wc -l|cut -d' ' -f1`
A_I_REDIportal=`cut -f 8,19 ${result}.editing|grep 'ed'|grep 'AG\|TC'|wc -l|cut -d' ' -f1`
echo -e "$raw\t$quality\t$rmsk\t$snp\t$splicesite\t$somatic\t$final\t$A_to_I\t$C_to_U\t$total_alu\t$total_rediportal\t$A_I_Alu\t$A_I_REDIportal" >>$out_path/${analysis_id}/${analysis_id}_filtering_summary.txt

