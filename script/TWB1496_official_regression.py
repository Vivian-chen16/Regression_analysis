#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv(r'/work1/viviane1695/tbb1496/updates/dataframe/tbb_updated_af_vqsr_cleaned.txt',sep='\t')

############################################################################### Filter distribution #####################################################################################################
df['tbbaf_all'] = df['tbbaf_all'].apply(pd.to_numeric, errors='coerce')
df['tbbaf_all'].fillna(0, inplace=True)

value = df['Filter'].value_counts(normalize=True) * 100
#print(np.round(value,2))

with open('/work1/viviane1695/tbb1496/updates/1496_official_filter_percent.txt', 'w') as f:
    for item in [np.round(value,2)]:
        f.write("%s\n" % item)

df['Filter'] = df['Filter'].map({'PASS': 'PASS(70.68%)', 
                                 'VQSRTrancheSNP99.80to99.90': 'SNV99.80-99.90(12.25%)',
                                'VQSRTrancheSNP99.70to99.80':'SNV99.70-99.80(9.90%)',
                                'VQSRTrancheSNP99.90to100.00':'SNV99.90-100.00(7.07%)',
                                'VQSRTrancheINDEL99.90to100.00':'INDEL99.90-100.00(0.08%)',
                                'VQSRTrancheINDEL99.70to99.80':'INDEL99.70-99.80(0.02%)'})

#### Pic 0 ####
fig, ax = plt.subplots()

sns.scatterplot(data=df, x='TWB1496', y='tbbaf_all', hue='Filter', hue_order= ['PASS(70.68%)', 'SNV99.90-100.00(7.07%)','SNV99.80-99.90(12.25%)','SNV99.70-99.80(9.90%)','INDEL99.90-100.00(0.08%)','INDEL99.70-99.80(0.02%)'], zorder=2,legend='full')
ax.legend(title='')

plt.xlabel('TWB1496', fontsize=12)
plt.ylabel('Taiwan official (tbbaf_all)', fontsize=12)
plt.grid(True)

plt.legend(title='VQSR Filter \n N=18,264', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0) #place legend outside top right corner of plot

#plt.show()
fig.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_official_filter.png', bbox_inches="tight", dpi=600)

############################################################################### population difference ###################################################################################################
df['Difference'] = np.where(df['TWB1496'] == df['tbbaf_all'], 0, df['TWB1496'] - df['tbbaf_all'])
df['Difference'] = df['Difference'].abs()

df['Filter_code'] = np.where(df['Filter']=='PASS(70.68%)',1, 2)

conditions = [
    (df['Difference']>0.05)&(df['Filter_code'] == 1),
    (df['Difference']>0.05)&(df['Filter_code'] == 2),
    (df['Difference']<=0.05)] # create a list of our conditions: Diff > 0.05

values = ['Sig_PASS','Sig_other', 'None'] # create a list to assign for each condition
df['Difference_Flag'] = np.select(conditions, values)

percent = df['Difference_Flag'].value_counts(normalize=True) * 100
#print(np.round(percent,2))

with open('/work1/viviane1695/tbb1496/updates/1496_official_DF_percent.txt', 'w') as f:
    for item in [np.round(percent,2)]:
        f.write("%s\n" % item)

df['Difference_Flag'] = df['Difference_Flag'].map({'None': 'Others(94.96%)', 
                                                   'Sig_other': 'Non-PASS(4.01%)',
                                                  'Sig_PASS':'PASS(1.03%)'})

#### Pic 1 ####
x = df['TWB1496']
y = df['tbbaf_all']

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
predict_y = intercept + slope * x

fig_1, ax_1 = plt.subplots()

sns.scatterplot(data=df, x='TWB1496', y='tbbaf_all', hue='Difference_Flag', hue_order= ['Others(94.96%)', 'Non-PASS(4.01%)','PASS(1.03%)'], zorder=2, palette=['silver','lightgreen','forestgreen'], legend='full')
ax_1.legend(title='')

plt.plot(x, predict_y, color='navy',label='$r^2 = {}$'.format(np.round(r_value**2,4)))
plt.plot(x, predict_y, color='navy',label='$y = {} + {}x$'.format(np.round(intercept,4), np.round(slope,4)))

plt.text(0,0.9,'$N=18,264$', fontsize=12) #print text on plot
plt.xlabel('TWB1496', fontsize=12)
plt.ylabel('Taiwan official(tbbaf_all)', fontsize=12)
plt.grid(True)

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0) #place legend outside top right corner of plot

#plt.show()
fig_1.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_official_DF_r2.png', bbox_inches="tight", dpi=600)

########################################################################### population difference of 'PASS' #############################################################################################

df1 = df[df['Filter']=='PASS(70.68%)']

df1['tbbaf_all'] = df1['tbbaf_all'].apply(pd.to_numeric, errors='coerce')
df1['tbbaf_all'].fillna(0, inplace=True)

df1['Difference'] = np.where(df1['TWB1496'] == df1['tbbaf_all'], 0, df1['TWB1496'] - df1['tbbaf_all'])
df1['Difference'] = df1['Difference'].abs()

conditions = [
    (df1['Difference']>0.05),
    (df1['Difference']<=0.05)] # create a list of our conditions: Diff > 0.05


values = ['Sig', '0'] # assign for each condition
df1['Difference_Flag'] = np.select(conditions, values)
percent1 = df1['Difference_Flag'].value_counts(normalize=True) * 100
#print(np.round(percent,2))

with open('/work1/viviane1695/tbb1496/updates/1496_official_DF_PASS.txt', 'w') as f:
    for item in [np.round(percent1,2)]:
        f.write("%s\n" % item)

df1['Difference_Flag'] = df1['Difference_Flag'].map({'0': 'Others(98.54%)', 'Sig': 'DF> 0.05 (1.46%)'})

#### Pic 3 ####
x1 = df1['TWB1496']
y1 = df1['tbbaf_all']

slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(x1,y1)
predict_y1 = intercept_1 + slope_1 * x1

fig_3, ax_3 = plt.subplots()

sns.scatterplot(data=df1, x='TWB1496', y='tbbaf_all', hue='Difference_Flag', hue_order= ['Others(98.54%)', 'DF> 0.05 (1.46%)'], zorder=2, palette=['silver','forestgreen'], legend='full')
ax_3.legend(title='')

plt.plot(x1, predict_y1, color='navy',label='$r^2 = {}$'.format(np.round(r_value_1**2,4)))
plt.plot(x1, predict_y1, color='navy',label='$y = {} + {}x$'.format(np.round(intercept_1,4), np.round(slope_1,4)))

plt.text(0,0.95,'$N=12,909$', fontsize=12) #print text on plot
plt.xlabel('TWB1496', fontsize=12)
plt.ylabel('Taiwan official(tbbaf_all)', fontsize=12)
plt.grid(True)

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)#place legend outside top right corner of plot

#plt.show()
fig_3.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_official_PASS_DF_r2.png', bbox_inches="tight", dpi=600)

########################################################################### population specific of 'PASS' #############################################################################################

conditions_eas = [
    (df1['TWB1496'] > 0.01) & (df1['tbbaf_all'] == 0),
    (df1['tbbaf_all'] > 0.01) & (df1['TWB1496'] == 0)]

values_eas = ['TWB1496', 'tbbaf_all']

df1['Specific_Flag'] = np.select(conditions_eas, values_eas)
specific = df1['Specific_Flag'].value_counts(normalize=True) * 100
#print(np.round(specific,2))

with open('/work1/viviane1695/tbb1496/updates/1496_official_specific_percent.txt', 'w') as f:
    for item in [np.round(specific,2)]:
        f.write("%s\n" % item)
#### Pic 4 ####
select_color = df1.loc[df1['Specific_Flag'] != '0']

x1 = df1['TWB1496']
y1 = df1['tbbaf_all']

slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(x1,y1)
predict_y1 = intercept_1 + slope_1 * x1

fig_4, ax_4 = plt.subplots()

sns.scatterplot(data=df1, x='TWB1496', y='tbbaf_all', hue='Specific_Flag', hue_order=['0','tbbaf_all','TWB1496'],palette=['silver','deeppink','dodgerblue'], zorder=1, legend=False)

select_color['Specific_Flag'] = select_color['Specific_Flag'].map({'tbbaf_all':'TWB official(0.02%)','TWB1496':'TWB1496(0.29%)'})
sns.scatterplot(data=select_color, x='TWB1496', y='tbbaf_all',hue='Specific_Flag', hue_order=['TWB official(0.02%)','TWB1496(0.29%)'], palette=['deeppink','dodgerblue'], zorder=2, legend='full')
ax_4.legend(title='')

plt.plot(x1, predict_y1, color='navy',label='$r^2 = {}$'.format(np.round(r_value_1**2,4)))
plt.plot(x1, predict_y1, color='navy',label='$y = {} + {}x$'.format(np.round(intercept_1,4), np.round(slope_1,4)))

plt.text(0,0.95,'$N=12,909$', fontsize=12) #print text on plot
plt.xlabel('TWB1496', fontsize=12)
plt.ylabel('Taiwan official(tbbaf_all)', fontsize=12)
plt.grid(True)

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0) #place legend outside top right corner of plot

#plt.show()
fig_4.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_official_PASS_specific_r2.png', bbox_inches="tight", dpi=600)


