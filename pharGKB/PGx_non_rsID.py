#!/usr/bin/env python3
# cofing:utf-8

import pandas as pd
import numpy as np

name = ['Chr', 'Start', 'Ref', 'Alt','Info','AF_Category']
df_AF = pd.read_csv(r'/work1/viviane1695/pharGKB/20210723/AF_category_extract_from_PGx_pos.txt', sep='\t',header=None, names=name, low_memory=False)


name1 = ['Chr', 'Start', 'End', 'Ref', 'variant']
df_pos = pd.read_csv(r'/work1/viviane1695/pharGKB/20210723/PGx_list_position_chr.txt', sep='\t',header=None, names=name1, low_memory=False)

df_phar = pd.read_csv(r'/work1/viviane1695/pharGKB/pharm_TWB_AF_20210609_Other.txt', sep='\t', low_memory=False)
df_other = df_phar[df_phar['variant'].str.contains('rs') == False]


df = df_pos.merge(df_phar)
df1 = df.merge(df_AF, how = 'left')
df1 = df1[['Chr','Start','End','Ref','Alt','variant','gene','type','level of evidence','chemicals','phenotypes','Info','AF_Category']]

df_non = df1[df1['Info'].isnull()]

df_final = pd.concat([df_non,df_other])
df_final = df_final.fillna(".")
df_final.to_csv(r'/work1/viviane1695/pharGKB/20210723/TWB_1496_PGx_table_non_rsID.txt', header=True, index=False, sep='\t')
