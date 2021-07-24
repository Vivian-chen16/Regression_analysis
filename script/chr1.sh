#!/bin/bash
#PBS -q ngs96G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N chr1_sub
#PBS -o /work1/viviane1695/tbb1496/updates/run.log/chr1_out
#PBS -e /work1/viviane1695/tbb1496/updates/run.log/chr1_err
#PBS -M vivianchen1695@gmail.com
#PBS -m a

awk -F "\t" '{if($1=="chr1") {print $0} else if($2=="Start") {print $0}}' OFS="\t" /project/GP1/j116831/1496_Annotation/anno_clinvar202105/TBB_1496_joing_calling.SNP_INDEL.recaled.decomposed.normalized.cut.hg19_multianno.txt | awk -v cols='Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc.refGene,AAChange.refGene,Func.knownGene,Gene.knownGene,GeneDetail.knownGene,ExonicFunc.knownGene,AAChange.knownGene,CLNALLELEID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG,avsnp150,tbbaf_all,tbbaf_illumina,gnomAD_2_genome_AF_eas,gnomAD_2_exome_AF_eas,tbbaf' 'BEGIN {
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
}' > /work1/viviane1695/tbb1496/updates/tbb_subset_chr1_20210618.txt
