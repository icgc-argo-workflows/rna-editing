#!/bin/bash
##############################################################################################
# Copyright(c) bjcancer, all rights reserved
# @file        : reditools_dnarna_filtering.sh
# @author      : yangyang
# @email       : yangyang@bjcancer.org
# @revision    : 2021-07-8 03:14:20
# @description : pipeline of reditools_dnarna annotation and filtering;
# @description : parament from reference: NATURE PROTOCOL (PMID:31996844) (DOI:10.1038/s41596-019-0279-7)
###############################################################################################
if [ $# != 10 ] ; then 
  echo -e "Usage: reditools_dnarna_filtering.sh reditools_path out_path analysis_id rna_bam dna_bam ref_fa rmsk_file snp151_file rediportal_file splicesites_file"
  echo -e "reditools_path: \n"
  echo -e "out_path: \n"
  echo -e "analysis_id: \n"
  echo -e "rna_bam: \n"
  echo -e "dna_bam: \n"
  echo -e "ref_fa: \n"
  echo -e "rmsk_file: \n"
  echo -e "snp151_file: \n"
  echo -e "rediportal_file: \n"
  echo -e "splicesites_file: \n"
exit 1;
fi

reditools_path=$1
out_path=$2
analysis_id=$3
rna_bam=$4
dna_bam=$5
ref_fa=$6
rmsk_file=$7
snp151_file=$8
rediportal_file=$9
splicesites_file=$10

### prepared files
# rmsk_file=/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/rmsk/rmsk.sorted.gtf.gz
# snp151_file=/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/snp151/snp151.sorted.gtf.gz
# rediportal_file=/mnt/dellfs/projects/RNA_editing/analysis_test/script/REDItools/resource/rediportal/atlas.gtf.gz
# splicesites_file=/mnt/dellfs/projects/RNA_editing/rna_editing_pipeline_sn/resources/Gencode_annotation/gencode.v30lift37.splicesites.txt

mkdir $out_path/${analysis_id}
cd $out_path/${analysis_id}
result=$out_path/${analysis_id}/outTable_${analysis_id}.out
threads=10

# 9. Enter the REDItoolDnaRna.py output folder:
in_file=`ls $reditools_path/DnaRna_${analysis_id}*/outTable_*`

# 10. Exclude invariant positions as well as positions not supported by ≥10 WGS reads:
## column8: ; column10: Pvalue ; 
awk 'FS="\t" {if ($8!="-" && $10>=10 && $13=="-") print}' $in_file >$result

# 11. Annotate positions using RepeatMasker and dbSNP annotations:
python /REDItools/accessory/AnnotateTable.py -a $rmsk_file -n rmsk -r chr -i $result -o ${result}.rmsk -u
python /REDItools/accessory/AnnotateTable.py -a $snp151_file -n snp151 -r chr -i ${result}.rmsk -o ${result}.rmsk.snp -u

# 12. Create a first set of positions selecting sites supported by at least four RNAseq reads and three mismatch:
python /REDItools/accessory/selectPositions.py -i ${result}.rmsk.snp -d -1 -c 4 -v 3 -f 0.1 -o ${result}.rmsk.snp.sel1

# 13. Create a second set of positions selecting sites supported by ≥10 RNAseq reads, three mismatches and minimum editing frequency of 0.1:
#python /REDItools/accessory/selectPositions.py -i ${result}.rmsk.snp -d -1 -c 10 -v 3 -f 0.1 -o ${result}.rmsk.snp.sel2

# 14. Select ALU sites from the first set of positions:
# awk 'FS="\t" {if ($1!="M" && substr($16,1,3)=="Alu" && $17=="-" && $8!="-" && $10>=10 && $13=="-") print}' ${result}.rmsk.snp.sel1 > ${result}.rmsk.snp.alu
awk 'FS="\t" {if ($1!="M" && substr($16,1,3)=="Alu" && $17=="-" && $8!="-" && $10>=10 && $13=="-") print}' ${result}.rmsk.snp.sel1 > ${result}.rmsk.snp.alu

# 15. Select REP NON ALU sites from the second set of positions, excluding sites in Simple repeats or Low complexity regions;
# ≥10 WGS reads and gFrequency<=0.05 and Frequency>=0.1
awk 'FS="\t" {if ($1!="M" && substr($16,1,3)!="Alu" && $15!="-" && $15!="Simple_repeat" && $15!="Low_complexity" && $17=="-" && $8!="-" && $10>=10 && $14<=0.05 && $9>=0.1) print}' ${result}.rmsk.snp.sel1 > ${result}.rmsk.snp.nonalu

# 16. Select NON REP sites from the second set of positions:
awk 'FS="\t" {if ($1!="M" && substr($16,1,3)!="Alu" && $15=="-" && $17=="-" && $8!="-" && $10>=10 && $14<=0.05 && $9>=0.1) print}' ${result}.rmsk.snp.sel1 > ${result}.rmsk.snp.nonrep

# 17. Annotate ALU, REP NON ALU and NON REP sites using known editing events from REDIportal:
python /REDItools/accessory/AnnotateTable.py -a $rediportal_file -r chr -n ed -k R -c 1 -i ${result}.rmsk.snp.alu -o ${result}.rmsk.snp.alu.ed -u
python /REDItools/accessory/AnnotateTable.py -a $rediportal_file -r chr -n ed -k R -c 1 -i ${result}.rmsk.snp.nonalu -o ${result}.rmsk.snp.nonalu.ed -u 
python /REDItools/accessory/AnnotateTable.py -a $rediportal_file -r chr -n ed -k R -c 1 -i ${result}.rmsk.snp.nonrep -o ${result}.rmsk.snp.nonrep.ed -u

# 18. Extract known editing events from ALU, REP NON ALU and NON REP sites:
mv ${result}.rmsk.snp.alu.ed $out_path/${analysis_id}/alu
mv ${result}.rmsk.snp.nonalu.ed $out_path/${analysis_id}/nonalu
mv ${result}.rmsk.snp.nonrep.ed $out_path/${analysis_id}/nonrep
cat $out_path/${analysis_id}/alu $out_path/${analysis_id}/nonalu $out_path/${analysis_id}/nonrep > $out_path/${analysis_id}/alu-nonalu-nonrep
awk 'FS="\t" {if ($19=="ed") print}' $out_path/${analysis_id}/alu-nonalu-nonrep > $out_path/${analysis_id}/knownEditing

# 19. Convert editing candidates in REP NON ALU and NON REP sites in GFF format for further filtering:
cat $out_path/${analysis_id}/nonalu $out_path/${analysis_id}/nonrep > $out_path/${analysis_id}/nonalu-nonrep
awk 'FS="\t" {if ($19!="ed") print}' $out_path/${analysis_id}/nonalu-nonrep > $out_path/${analysis_id}/pos.txt
# where $19!="ed" selects novel RNA editing events.
python /REDItools/accessory/TableToGFF.py -i $out_path/${analysis_id}/pos.txt -s -t -o $out_path/${analysis_id}/pos.gff

# 20 Convert editing candidates in ALU sites in GFF format for further filtering:
awk 'FS="\t" {if ($19!="ed") print}' $out_path/${analysis_id}/alu > $out_path/${analysis_id}/posalu.txt
python /REDItools/accessory/TableToGFF.py -i $out_path/${analysis_id}/posalu.txt -s -t -o $out_path/${analysis_id}/posalu.gff

# 21. Launch REDItoolDnaRna.py on ALU sites using stringent criteria to recover potential editing candidates:
python /REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t $threads -i $rna_bam -f $ref_fa \
-c 10,4 -q 30,30 -m 30,30 -p -u -a 11-6 -l -v 3 -n 0.1 \
-e -T $out_path/${analysis_id}/posalu.sorted.gff.gz -w $splicesites_file -R -o $out_path/${analysis_id}/firstalu

# 22. Launch REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria to recover RNAseq reads harboring reference mismatches:
python /REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t $threads -i $rna_bam -f $ref_fa \
-c 10,4 -q 30,30 -m 30,30 -p -u -a 11-6 -l -v 3 -n 0.1 \
-e -T $out_path/${analysis_id}/pos.sorted.gff.gz -w $splicesites_file --reads -R --addP -o $out_path/${analysis_id}/first

# 23. Launch pblat on RNAseq reads harboring reference mismatches from Step 22 and select multi-mapping reads:
file_outReads=`ls $out_path/${analysis_id}/first/DnaRna_*/outReads_*`
/mnt/dellfs/pub/tools/pblat/icebert-pblat-98df9e0/pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 \
-minIdentity=0 $ref_fa $file_outReads $out_path/${analysis_id}/reads.psl
python /REDItools/accessory/readPsl.py $out_path/${analysis_id}/reads.psl $out_path/${analysis_id}/badreads.txt

# 24. Extract RNAseq reads harboring reference mismatches from Step 22 and remove duplicates:
cd $out_path/${analysis_id}
first_outPosReads=`ls $out_path/${analysis_id}/first/DnaRna_*/outPosReads_*`
sort -k1,1 -k2,2n -k3,3n $first_outPosReads | mergeBed > $out_path/${analysis_id}/bed
samtools view -@ 4 -L bed -h -b $rna_bam > $out_path/${analysis_id}/${analysis_id}_bed.bam
samtools sort -@ 4 -n $out_path/${analysis_id}/${analysis_id}_bed.bam -o $out_path/${analysis_id}/${analysis_id}_bed_ns.bam
samtools fixmate -@ 4 -m $out_path/${analysis_id}/${analysis_id}_bed_ns.bam $out_path/${analysis_id}/${analysis_id}_bed_ns_fx.bam
samtools sort -@ 4 $out_path/${analysis_id}/${analysis_id}_bed_ns_fx.bam -o $out_path/${analysis_id}/${analysis_id}_bed_ns_fx_st.bam
samtools markdup -r -@ 4 $out_path/${analysis_id}/${analysis_id}_bed_ns_fx_st.bam $out_path/${analysis_id}/${analysis_id}_bed_dedup.bam
samtools index $out_path/${analysis_id}/${analysis_id}_bed_dedup.bam

# 25. Re-run REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria, deduplicated reads and mis-mapping info:
python /REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t $threads -i $out_path/${analysis_id}/${analysis_id}_bed_dedup.bam \
-f $ref_fa -c 10,4 -q 30,30 -m 30,30 -p -u -a 11-6 -l -v 3 -n 0.1 \
-e -T $out_path/${analysis_id}/pos.sorted.gff.gz -w $splicesites_file -R -b $out_path/${analysis_id}/badreads.txt --rmIndels -o $out_path/${analysis_id}/second

# 26. Collect filtered ALU, REP NON ALU and NON REP sites:
cd $out_path/${analysis_id}
python /REDItools/NPscripts/collect_editing_candidates.py
sort -k1,1 -k2,2n $out_path/${analysis_id}/editing.txt > $out_path/${analysis_id}/editing_sorted.txt

# 27. Inspect the distribution of editing candidates to look at A-to-I enrichment:
python /REDItools/NPscripts/get_Statistics.py

# summary of editing_sorted.txt
echo -e  "total\tA_to_I\tC_to_U\ttotal_Alu\ttotal_REDIportal\tA_I_Alu\tA_I_REDIportal" >$out_path/${analysis_id}/${analysis_id}_filtered_hc_anno_summary.txt
total=`wc -l $out_path/${analysis_id}/editing_sorted.txt|cut -d' ' -f1`
A_to_I=`grep 'AG\|TC' $out_path/${analysis_id}/editing_sorted.txt|wc -l|cut -d' ' -f1`
C_to_U=`grep 'GA\|CT' $out_path/${analysis_id}/editing_sorted.txt|wc -l|cut -d' ' -f1`
total_Alu=`grep Alu $out_path/${analysis_id}/editing_sorted.txt|wc -l|cut -d' ' -f1`
total_noAlu=`grep -v Alu $out_path/${analysis_id}/editing_sorted.txt|wc -l|cut -d' ' -f1`
total_REDIportal=`grep -w ed $out_path/${analysis_id}/editing_sorted.txt|wc -l|cut -d' ' -f1`
A_I_Alu=`grep Alu $out_path/${analysis_id}/editing_sorted.txt|grep 'AG\|TC'|wc -l|cut -d' ' -f1`
A_I_REDIportal=`grep -w ed $out_path/${analysis_id}/editing_sorted.txt|grep 'AG\|TC'|wc -l|cut -d' ' -f1`
echo -e "$total\t$A_to_I\t$C_to_U\t$total_Alu\t$total_REDIportal\t$A_I_Alu\t$A_I_REDIportal" >>$out_path/${analysis_id}/${analysis_id}_filtered_hc_anno_summary.txt

# summary of alu-nonalu-nonrep
echo -e  "total\tA_to_I\tC_to_U\ttotal_Alu\ttotal_REDIportal\tA_I_Alu\tA_I_REDIportal" >$out_path/${analysis_id}/${analysis_id}_filtered_anno_summary.txt
total=`wc -l $out_path/${analysis_id}/alu-nonalu-nonrep|cut -d' ' -f1`
A_to_I=`grep 'AG\|TC' $out_path/${analysis_id}/alu-nonalu-nonrep|wc -l|cut -d' ' -f1`
C_to_U=`grep 'GA\|CT' $out_path/${analysis_id}/alu-nonalu-nonrep|wc -l|cut -d' ' -f1`
total_Alu=`grep Alu $out_path/${analysis_id}/alu-nonalu-nonrep|wc -l|cut -d' ' -f1`
total_noAlu=`grep -v Alu $out_path/${analysis_id}/alu-nonalu-nonrep|wc -l|cut -d' ' -f1`
total_REDIportal=`grep -w ed $out_path/${analysis_id}/alu-nonalu-nonrep|wc -l|cut -d' ' -f1`
A_I_Alu=`grep Alu $out_path/${analysis_id}/alu-nonalu-nonrep|grep 'AG\|TC'|wc -l|cut -d' ' -f1`
A_I_REDIportal=`grep -w ed $out_path/${analysis_id}/alu-nonalu-nonrep|grep 'AG\|TC'|wc -l|cut -d' ' -f1`
echo -e "$total\t$A_to_I\t$C_to_U\t$total_Alu\t$total_REDIportal\t$A_I_Alu\t$A_I_REDIportal" >>$out_path/${analysis_id}/${analysis_id}_filtered_anno_summary.txt

