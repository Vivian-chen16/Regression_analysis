#!/bin/bash
#PBS -q ngs96G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N tbb_condition
#PBS -o /work1/viviane1695/tbb1496/updated_annotation/run.log/condition_20210618_out
#PBS -e /work1/viviane1695/tbb1496/updated_annotation/run.log/condition_20210618_err
#PBS -M vivianchen1695@gmail.com
#PBS -m a

awk '{ if($11=='ebonic')
{
 print $0
}
else if($11 ~ /exonic/)
{
 print $0
}
else if($11 ~ /ncRNA_exonic/)
{
 print $0
}
else if($11=='splicing')
{
 print $0
}
else if($11=='ncRNA_splicing')
{
 print $0
}
}' /work1/viviane1695/tbb1496/updated_annotation/tbb_subset_20210618.txt > /work1/viviane1695/tbb1496/updated_annotation/tbb_cond1_20210618.txt

awk '{ if($14 ~/nonsynonymous/)
{
 print $0
}
else if($14 ~/frameshift/)
{
 print $0
}
else if($14=='stopgain')
{
 print $0
}
else if($14 ~/nonframeshift/)
{
 print $0
}
else if($14=='stoploss') 
{
 print $0
}
}' /work1/viviane1695/tbb1496/updated_annotation/tbb_cond1_20210618.txt > /work1/viviane1695/tbb1496/updated_annotation/tbb_cond2_20210618.txt

echo -e 'Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tFunc.knownGene\tGene.knownGene\tGeneDetail.knownGene\tExonicFunc.knownGene\tAAChange.knownGene\tCLNALLELEID\tCLNDN\tCLNDISDB\tCLNREVSTAT\tCLNSIG\tavsnp150\ttbbaf_all\ttbbaf_illumina\tgnomAD_2_genome_AF_eas\tgnomAD_2_exome_AF_eas\ttbbaf' | cat - /work1/viviane1695/tbb1496/updated_annotation/tbb_cond2_20210618.txt > /work1/viviane1695/tbb1496/updated_annotation/tbb_subset_conditional_20210618.txt

