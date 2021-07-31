#!/usr/bin/env python3
# cofing:utf-8

import pandas as pd
import dask.dataframe as dd
import numpy as np

df_1496 = pd.read_csv(r'/work1/viviane1695/pharGKB/20210723/tbb_subset_20210618_without_tbbaf.txt', sep='\t', low_memory=False)
df_1496 = df_1496.rename(columns={'avsnp150':'variant'})

name = ['Chr', 'Start', 'Ref', 'Alt','Info','AF_Category']
df_AF = pd.read_csv(r'/work1/viviane1695/pharGKB/20210723/AF_category_extract_from_PGx_pos.txt', sep='\t',header=None, names=name, low_memory=False)


name1 = ['Chr', 'Start', 'End', 'Ref', 'variant']
df_pos = pd.read_csv(r'/work1/viviane1695/pharGKB/20210723/PGx_list_position_chr.txt', sep='\t',header=None, names=name1, low_memory=False)
df_pos.drop_duplicates(keep=False,inplace=True)

columns = ['Chr', 'Start', 'ID', 'Ref', 'Alt', 'Quality', 'Filter']
df_VQSR = pd.read_csv(r'/work1/viviane1695/pharGKB/20210723/VQSR_PGx_list.txt', header=None, names=columns, sep='\t', low_memory=False)
df_VQSR = df_VQSR.drop(columns=['ID'])

df_phar = pd.read_csv(r'/work1/viviane1695/pharGKB/pharm_TWB_AF_20210609_Other.txt', sep='\t', low_memory=False)


df = df_pos.merge(df_phar)
df1 = df.merge(df_AF, how = 'left')
df2 = df1.merge(df_VQSR, how = 'left')
df2 = df2[['Chr','Start','End','Ref','Alt','variant','gene','type','level of evidence','chemicals','phenotypes','Quality','Filter','Info','AF_Category']]

df_final = df2.merge(df_1496, how='left')
df_final = df_final[df_final['Info'].isnull() == False]

df_final.to_csv(r'/work1/viviane1695/pharGKB/20210723/TWB_1496_updated_annotation_final_PGx_table.txt', header=True, index=False, sep='\t')
