#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.lines import Line2D

df = pd.read_csv(r'/work1/viviane1695/tbb1496/updates/dataframe/tbb_updated_af_vqsr_cleaned.txt',sep='\t')

df['Filter_code'] = df['Filter'].str.split('to').str[1]
df['Filter_code'].fillna('PASS',inplace=True)

df998 = df[df['Filter_code']=='99.80']
df999 = df[df['Filter_code']=='99.90']

df1 = df[df['Filter_code']=='PASS']
df2 = pd.concat([df1, df998], ignore_index=True)
df3 = pd.concat([df2, df999], ignore_index=True)

value_1 = df['Filter_code'].value_counts(normalize=True) * 100
#print(np.round(value_1,2))

with open('/work1/viviane1695/tbb1496/updates/from_4plot.txt', 'w') as f:
    for item in [np.round(value_1,2)]:
        f.write("%s\n" % item)

df['Filter_code'] = df['Filter_code'].map({'PASS': 'PASS(70.68%)', '99.80': '99.70-99.80(9.92%)','99.90':'99.80-99.90(12.25%)','100.00':'99.90-100.00(7.15%)'})
df1['Filter_code'] = df1['Filter_code'].map({'PASS': 'PASS(70.68%)', '99.80': '99.70-99.80(9.92%)','99.90':'99.80-99.90(12.25%)','100.00':'99.90-100.00(7.15%)'})
df2['Filter_code'] = df2['Filter_code'].map({'PASS': 'PASS(70.68%)', '99.80': '99.70-99.80(9.92%)','99.90':'99.80-99.90(12.25%)','100.00':'99.90-100.00(7.15%)'})
df3['Filter_code'] = df3['Filter_code'].map({'PASS': 'PASS(70.68%)', '99.80': '99.70-99.80(9.92%)','99.90':'99.80-99.90(12.25%)','100.00':'99.90-100.00(7.15%)'})

#For TWB1496
x = df['TWB1496']  
y = df['gnomAD_2_Final']
  
x1 = df1['TWB1496']
y1 = df1['gnomAD_2_Final']
x2 = df2['TWB1496']
y2 = df2['gnomAD_2_Final']
x3 = df3['TWB1496']
y3 = df3['gnomAD_2_Final']

colors = {'PASS(70.68%)':'tab:blue', '99.70-99.80(9.92%)':'tab:green', '99.80-99.90(12.25%)':'tab:orange', '99.90-100.00(7.15%)':'tab:red'}

#Initialise the subplot function using number of rows and columns
figure, axis = plt.subplots(2, 2)
  
#For PASS
slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(x1,y1)
predict_y1 = intercept_1 + slope_1 * x1
label_1 = '$r^2 = {}$'.format(np.round(r_value_1**2,4))

axis[0,0].scatter(x1, y1, c = df1['Filter_code'].map(colors), s=12, edgecolors='w', linewidth=0.3) #c=lable color; s=dot size
axis[0,0].plot(x1, predict_y1, color='navy')
axis[0,0].text(0,0.9,label_1, fontsize=10)#print text on plot
axis[0,0].grid(True)

#For PASS + 99.70-99.80
slope_2, intercept_2, r_value_2, p_value_2, std_err_2 = stats.linregress(x2,y2)
predict_y2 = intercept_2 + slope_2 * x2
label_2 = '$r^2 = {}$'.format(np.round(r_value_2**2,4))

axis[0,1].scatter(x2, y2, c = df2['Filter_code'].map(colors),s=12, edgecolors='w', linewidth=0.3)
axis[0,1].plot(x2, predict_y2, color='navy')
axis[0,1].text(0,0.9,label_2, fontsize=10) #print text on plot
axis[0,1].grid(True)
handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=6) for k, v in colors.items()]
axis[0,1].legend(title='VQSR Filter \n N=18,624', handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')

#For PASS + 99.70-99.80 + 99.80-99.90
slope_3, intercept_3, r_value_3, p_value_3, std_err_3 = stats.linregress(x3,y3)
predict_y3 = intercept_3 + slope_3 * x3
label_3 = '$r^2 = {}$'.format(np.round(r_value_3**2,4))

axis[1,0].scatter(x3, y3, c = df3['Filter_code'].map(colors), s=12, edgecolors='w', linewidth=0.3) #c=lable color; s=dot size
axis[1,0].plot(x3, predict_y3, color='navy')
axis[1,0].text(0,0.9,label_3, fontsize=10)#print text on plot
axis[1,0].grid(True)

#For ALL
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
predict_y = intercept + slope * x
label = '$r^2 = {}$'.format(np.round(r_value**2,4))

axis[1,1].scatter(x, y, c = df['Filter_code'].map(colors), s=12, edgecolors='w',linewidth=0.3) #c=lable color; s=dot size
axis[1,1].plot(x, predict_y, color='navy')
axis[1,1].text(0,0.9,label, fontsize=10)#print text on plot
axis[1,1].grid(True)

#figure.suptitle('East Asian vs TWB1496, MAF in different quality level')
figure.text(0.5, 0.0, 'TWB1496', ha='center', fontsize=12)
figure.text(0.0, 0.5, 'East Asian(EAS, gnomAD)', va='center', rotation='vertical', fontsize=12)
#plt.show()

figure.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_gnomAD_fourplot_r2.png', bbox_inches="tight", dpi=600)
