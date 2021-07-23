#!/bin/bash
#PBS -q ngs192G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N Merge
#PBS -o /work1/viviane1695/tbb1496/updated_annotation/run.log/merge_out
#PBS -e /work1/viviane1695/tbb1496/updated_annotation/run.log/merge_err
#PBS -M vivianchen1695@gmail.com
#PBS -m a

export PATH=/pkg/biology/Python/Python3_default/bin:$PATH
cd /work1/viviane1695/tbb1496/updated_annotation

awk -F "\t" 'FNR == NR {strings[$1,$2]} NR > FNR && (($1,$2) in strings)' /work1/viviane1695/tbb1496/updated_annotation/tbb_subset_conditional_20210618.txt /project/GP1/j116831/1496_Annotation/RefAllele_MAF/20210515_interval_test/chrM_1_22_XY_MAF/20210611_dropGT/AF_category_final.txt | awk -F "\t" '{if($9 ~/A/) {print $1,$2,$4,$5,$9}}' OFS="\t" > /work1/viviane1695/tbb1496/updated_annotation/AF_Category_subset_nonsynonym.txt

python3 /work1/viviane1695/tbb1496/updated_annotation/tbb_af_vqsr_merge_and_clean.py

