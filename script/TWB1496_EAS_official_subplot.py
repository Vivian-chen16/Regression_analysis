#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.lines import Line2D

df = pd.read_csv(r'/work1/viviane1695/tbb1496/updated_annotation/dataframe/tbb_updated_af_vqsr_cleaned.txt',sep='\t')

df['tbbaf_all'] = df['tbbaf_all'].apply(pd.to_numeric, errors='coerce')
df['tbbaf_all'].fillna(0, inplace=True)

value = df['Filter'].value_counts(normalize=True) * 100
print(np.round(value,2))

df['Filter'] = df['Filter'].map({'PASS': 'PASS(70.68%)', 
                                 'VQSRTrancheSNP99.80to99.90': 'SNV99.80-99.90(12.25%)',
                                'VQSRTrancheSNP99.70to99.80':'SNV99.70-99.80(9.90%)',
                                'VQSRTrancheSNP99.90to100.00':'SNV99.90-100.00(7.07%)',
                                'VQSRTrancheINDEL99.90to100.00':'INDEL99.90-100.00(0.08%)',
                                'VQSRTrancheINDEL99.70to99.80':'INDEL99.70-99.80(0.02%)'})


y = df['gnomAD_2_Final']
  
x1 = df['TWB1496']
x2 = df['tbbaf_all']
colors = {'PASS(70.68%)':'tab:blue', 'SNV99.70-99.80(9.90%)':'tab:orange', 'SNV99.80-99.90(12.25%)':'tab:green', 'SNV99.90-100.00(7.07%)':'tab:red', 'INDEL99.70-99.80(0.02%)':'tab:purple', 'INDEL99.90-100.00(0.08%)':'tab:pink'}

# Initialise the subplot function using number of rows and columns
figure, axis = plt.subplots(1, 2)
  
# For TWB1496
slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(x1,y)
predict_y1 = intercept_1 + slope_1 * x1
label_10 = '$r^2 = {}$'.format(np.round(r_value_1**2,4))
label_11 = '$y = {} + {}x$'.format(np.round(intercept_1,4), np.round(slope_1,4))

axis[0].scatter(x1, y, c = df['Filter'].map(colors), s=20, edgecolors='w',linewidth=0.5) #c=lable color; s=dot size
axis[0].plot(x1, predict_y1, color='navy')
axis[0].set_xlabel('TWB1496', fontsize=12)
axis[0].set_ylabel('East Asian(EAS, gnomAD)', fontsize=12)
axis[0].text(0,1.1,label_10, fontsize=10)#print text on plot
axis[0].text(0,1.18,label_11, fontsize=10)
axis[0].grid(True)

  
# For Official
slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = stats.linregress(x2,y)
predict_y2 = intercept_2 + slope_2 * x2
label_21 = '$r^2 = {}$'.format(np.round(r_value_2**2,4))
label_22 = '$y = {} + {}x$'.format(np.round(intercept_2,4), np.round(slope_2,4))

axis[1].scatter(x2, y, c = df['Filter'].map(colors),s=20, edgecolors='w',linewidth=0.5)
axis[1].plot(x2, predict_y2, color='navy')
axis[1].set_xlabel('TWB Official(tbbaf_all)', fontsize=12)
axis[1].text(0,1.1,label_21, fontsize=10) #print text on plot
axis[1].text(0,1.18,label_22, fontsize=10)
axis[1].grid(True)
handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=6) for k, v in colors.items()]
axis[1].legend(title='VQSR Filter \n N=18,264', handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')


#plt.show()
figure.savefig('/work1/viviane1695/tbb1496/updated_annotation/figure/eas_twb1496_official_fillna_r2.png', bbox_inches="tight", dpi=600)

