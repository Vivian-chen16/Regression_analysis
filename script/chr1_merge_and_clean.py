#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import dask.dataframe as dd
import numpy as np

df_twb = pd.read_csv(r'/work1/viviane1695/tbb1496/updates/tbb_subset_chr1_20210618.txt', sep='\t', low_memory=False)

df_AF = pd.read_csv(r'/work1/viviane1695/tbb1496/updates/AF_Category_subset_chr1.txt', sep='\t', low_memory=False)
df_AF.columns = ['Chr','Start','Ref','Alt','Info','AF/MAF_Category']
#df_AF = df_AF.astype({'Start':'int64','Quality':'float64'})


df_VQSR = pd.read_csv(r'/project/GP1/j116831/1496_Annotation/RefAllele_MAF/20210515_interval_test/chrM_1_22_XY_MAF/20210611_dropGT/AF_VQSR_extraction.txt', sep='\t', low_memory=False)
df_VQSR.columns = ['Chr','Start','ID','Ref','Alt','Quality','Filter']
#df_VQSR = df_VQSR.astype({'Start':'int64'})

df1 = pd.merge(df_twb, df_VQSR)
df2 = pd.merge(df1, df_AF)
df2.to_csv(r'/work1/viviane1695/tbb1496/updates/dataframe/chr1_updated_af_vqsr_merge.txt', index=False, header=True, sep='\t', mode='a')

###################################################### Clean merged dataframe for regression/scatter plot ########################################################################################

df = pd.read_csv(r'/work1/viviane1695/tbb1496/updates/dataframe/chr1_updated_af_vqsr_merge.txt', sep='\t', low_memory=False)
df['TWB1496'] = df['AF/MAF_Category'].str.split('/').str[0]

#calcalate final gnomAD AF
df[['TWB1496','gnomAD_2_genome_AF_eas', 'gnomAD_2_exome_AF_eas']] = df[['TWB1496','gnomAD_2_genome_AF_eas', 'gnomAD_2_exome_AF_eas']].apply(pd.to_numeric, errors='coerce')
df = df.dropna(subset=['gnomAD_2_genome_AF_eas', 'gnomAD_2_exome_AF_eas'], how='all')
df.reset_index(drop=True, inplace=True)

def take_bigger(row):
    if row['gnomAD_2_genome_AF_eas'] == 'NaN':
        return row['gnomAD_2_exome_AF_eas'] 
    elif ((row['gnomAD_2_genome_AF_eas']!='NaN')&(row['gnomAD_2_exome_AF_eas']!='NaN')):
        if row['gnomAD_2_genome_AF_eas'] > row['gnomAD_2_exome_AF_eas']:
            return row['gnomAD_2_genome_AF_eas']
        else:
            return row['gnomAD_2_exome_AF_eas']
    else:
        return row['gnomAD_2_genome_AF_eas']

df['gnomAD_2_Final'] = df.apply (lambda row: take_bigger(row), axis=1)

df_1 = df[df['gnomAD_2_Final'].notna()]
df_2 = df[df['gnomAD_2_Final'].isna()]
df_2['gnomAD_2_Final'] = df_2['gnomAD_2_genome_AF_eas'].values

df = pd.concat([df_1,df_2], ignore_index=True)

#cleand data for regression analysis
df = df.sort_values(by=['Chr','Start'])
df.to_csv(r'/work1/viviane1695/tbb1496/updates/dataframe/chr1_updated_af_vqsr_cleaned.txt', index=False, header=True, sep='\t', mode='a')
