#!/bin/bash
#PBS -q ngs192G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N chr1_regression
#PBS -o /work1/viviane1695/tbb1496/updated_annotation/run.log/chr1_regression_out
#PBS -e /work1/viviane1695/tbb1496/updated_annotation/run.log/chr1_regression_err
#PBS -M vivianchen1695@gmail.com
#PBS -m a

export PATH=/pkg/biology/Python/Python3_default/bin:$PATH
cd /work1/viviane1695/tbb1496/updated_annotation

awk -F "\t" '{if($1 == "chr1") {print $1,$2,$4,$5,$8,$9}}' OFS="\t" /project/GP1/j116831/1496_Annotation/RefAllele_MAF/20210515_interval_test/chrM_1_22_XY_MAF/20210611_dropGT/AF_category_final.txt | awk -F "\t" '{if($6 ~/A/) {print $0}}' > /work1/viviane1695/tbb1496/updated_annotation/AF_Category_subset_chr1.txt

python3 chr1_merge_and_clean.py
python3 4subplot_by_filter_chr1.py
