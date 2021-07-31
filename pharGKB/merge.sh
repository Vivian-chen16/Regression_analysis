#!/bin/bash
#PBS -q ngs384G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N PGx_subset
#PBS -o /work1/viviane1695/pharGKB/run.log/merge0723_out
#PBS -e /work1/viviane1695/pharGKB/run.log/merge0723_err
#PBS -M vivianchen1695@gmail.com
#PBS -m a


# rsID to position 
awk -F "\t" 'FNR == NR {strings[$1]} NR > FNR && ($6 in strings)' OFS="\t" /work1/viviane1695/pharGKB/pharm_TWB_AF_20210609_Other.txt /project/GP1/u3710062/AI_SHARE/reference/annovar_2016Feb01/humandb/hg19_avsnp150.txt > /work1/viviane1695/pharGKB/20210723/PGx_list_position.txt

# subset needed column and sync values $1
awk -F "\t" '{print $1,$2,$3,$4,$6}' OFS="\t" /work1/viviane1695/pharGKB/20210723/PGx_list_position.txt | awk -F "\t" '$1="chr"$1' OFS="\t" > /work1/viviane1695/pharGKB/PGx_list_position_chr.txt

# position to af
awk -F "\t" 'FNR == NR {strings[$1,$2]} NR > FNR && (($1,$2) in strings)' /work1/viviane1695/pharGKB/20210723/PGx_list_position_chr.txt /project/GP1/j116831/1496_Annotation/RefAllele_MAF/20210515_interval_test/chrM_1_22_XY_MAF/20210611_dropGT/AF_category_final.txt |awk -F "\t" '{print $1,$2,$4,$5,$8,$9}' OFS="\t" > /work1/viviane1695/pharGKB/20210723/AF_category_extract_from_PGx_pos.txt

# extract VQSR for Ref heter
awk -F "\t" 'FNR == NR {strings[$1,$2]} NR > FNR && (($1,$2) in strings)' OFS="\t" /work1/viviane1695/pharGKB/20210723/PGx_list_position_chr.txt /project/GP1/j116831/1496_Annotation/RefAllele_MAF/20210515_interval_test/chrM_1_22_XY_MAF/20210611_dropGT/AF_VQSR_extraction.txt > /work1/viviane1695/pharGKB/20210723/VQSR_PGx_list.txt



######### updated: (1) remove twb1497 official data (2) separate into 2 part: PGx list within/without twb1496  ####################

awk -v cols='Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc.refGene,AAChange.refGene,Func.knownGene,Gene.knownGene,GeneDetail.knownGene,ExonicFunc.knownGene,AAChange.knownGene,CLNALLELEID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG,avsnp150,gnomAD_2_genome_AF_eas,gnomAD_2_exome_AF_eas' 'BEGIN {
   FS=OFS="\t"
   nc=split(cols, a, ",")
}
NR==1 {
   for (i=1; i<=NF; i++)
      hdr[$i]=i
}
{
   for (i=1; i<=nc; i++)
      if (a[i] in hdr)
         printf "%s%s", $hdr[a[i]], (i<nc?OFS:ORS)
}' /work1/viviane1695/tbb1496/updated_annotation/tbb_subset_20210618.txt > /work1/viviane1695/pharGKB/20210723/tbb_subset_20210618_without_tbbaf.txt

export PATH=/pkg/biology/Python/Python3_default/bin:$PATH
cd /work1/viviane1695/pharGKB/20210723

python3 PGx_merge_0723.py
python3 PGx_non_rsID.py
