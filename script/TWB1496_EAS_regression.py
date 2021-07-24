#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv(r'/work1/viviane1695/tbb1496/updates/dataframe/tbb_updated_af_vqsr_cleaned.txt',sep='\t')

###################################################################################### Filter distribution #############################################################################################

value = df['Filter'].value_counts(normalize=True) * 100
#print(np.round(value,2))

with open('/work1/viviane1695/tbb1496/updates/1496_EAS_filter_percent.txt', 'w') as f:
    for item in [np.round(value,2)]:
        f.write("%s\n" % item)

df['Filter'] = df['Filter'].map({'PASS': 'PASS(70.68%)', 
                                 'VQSRTrancheSNP99.80to99.90': 'SNV99.80-99.90(12.25%)',
                                'VQSRTrancheSNP99.70to99.80':'SNV99.70-99.80(9.90%)',
                                'VQSRTrancheSNP99.90to100.00':'SNV99.90-100.00(7.07%)',
                                'VQSRTrancheINDEL99.90to100.00':'INDEL99.90-100.00(0.08%)',
                                'VQSRTrancheINDEL99.70to99.80':'INDEL99.70-99.80(0.02%)'})
#### Pic 0 ####
fig_0, ax_0 = plt.subplots()
sns.scatterplot(data=df, x='TWB1496', y='gnomAD_2_Final', hue='Filter', hue_order= ['PASS(70.68%)', 'SNV99.90-100.00(7.07%)','SNV99.80-99.90(12.25%)','SNV99.70-99.80(9.90%)','INDEL99.90-100.00(0.08%)','INDEL99.70-99.80(0.02%)'], zorder=2,legend='full')
ax_0.legend(title='')

plt.xlabel('TWB1496', fontsize=12)
plt.ylabel('East Asian(EAS, gnomAD)', fontsize=12)
plt.grid(True)

plt.legend(title='VQSR Filter \n N=18,264', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

#plt.show()
fig_0.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_gnomAD_filter.png', bbox_inches="tight", dpi=600)


###################################################################################### population difference ###########################################################################################

df['Difference'] = np.where(df['TWB1496'] == df['gnomAD_2_Final'], 0, df['TWB1496'] - df['gnomAD_2_Final'])
df['Difference'] = df['Difference'].abs()

df['Filter_code'] = np.where(df['Filter']=='PASS(70.68%)',1, 2)
# create a list of conditions: Diff > 0.05
conditions = [
    (df['Difference']>0.05)&(df['Filter_code'] == 1),
    (df['Difference']>0.05)&(df['Filter_code'] == 2),
    (df['Difference']<=0.05)]

# create a list to assign for each condition
values = ['Sig_PASS','Sig_other', 'None']
df['Difference_Flag'] = np.select(conditions, values)
percent = df['Difference_Flag'].value_counts(normalize=True) * 100
#print(np.round(percent,2))

with open('/work1/viviane1695/tbb1496/updates/1496_EAS_DF_percent.txt', 'w') as f:
    for item in [np.round(percent,2)]:
        f.write("%s\n" % item)

df['Difference_Flag'] = df['Difference_Flag'].map({'None': 'Others(96.68%)', 
                                                   'Sig_other': 'Non-PASS(2.97%)',
                                                  'Sig_PASS':'PASS(0.34%)'})

#### Pic 1 ####
x = df['TWB1496']
y = df['gnomAD_2_Final']

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
predict_y = intercept + slope * x

fig_1, ax_1 = plt.subplots()

sns.scatterplot(data=df, x='TWB1496', y='gnomAD_2_Final', hue='Difference_Flag', hue_order= ['Others(96.68%)', 'Non-PASS(2.97%)','PASS(0.34%)'], zorder=2, palette=['silver','salmon','firebrick'], legend='full')
ax_1.legend(title='')

plt.plot(x, predict_y, color='navy',label='$r^2 = {}$'.format(np.round(r_value**2,4)))
plt.plot(x, predict_y, color='navy',label='$y = {} + {}x$'.format(np.round(intercept,4), np.round(slope,4)))

plt.text(0,0.9,'$N=18,624$', fontsize=12) #print text on plot
plt.xlabel('TWB1496', fontsize=12)
plt.ylabel('East Asian(EAS, gnomAD)', fontsize=12)
plt.grid(True)

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

#plt.show()
fig_1.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_gnomAD_DF_r2.png', bbox_inches="tight", dpi=600)

############################################################################### population difference for 'PASS' ########################################################################################
# ALL difference are non-pass 
# VQSR=='PASS
df1 = df[df['Filter']=='PASS(70.68%)']

df1['Difference'] = np.where(df1['TWB1496'] == df1['gnomAD_2_Final'], 0, df1['TWB1496'] - df1['gnomAD_2_Final'])
df1['Difference'] = df1['Difference'].abs()


# create a list of our conditions: Diff > 0.05
conditions = [
    (df1['Difference']>0.05),
    (df1['Difference']<=0.05),]

# create a list to assign for each condition
values = ['Sig', '0']
df1['Difference_Flag'] = np.select(conditions, values)
percent1 = df1['Difference_Flag'].value_counts(normalize=True) * 100
#print(np.round(percent1,2))

with open('/work1/viviane1695/tbb1496/updates/1496_EAS_DF_PASS.txt', 'w') as f:
    for item in [np.round(percent1,2)]:
        f.write("%s\n" % item)

df1['Difference_Flag'] = df1['Difference_Flag'].map({'0': 'Others(99.51%)', 'Sig': 'DF > 0.05 (0.49%)'})

#### Pic 2 ####
x1 = df1['TWB1496']
y1 = df1['gnomAD_2_Final']

slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(x1,y1)
predict_y1 = intercept1 + slope1 * x1

fig_2, ax_2 = plt.subplots()

sns.scatterplot(data=df1, x='TWB1496', y='gnomAD_2_Final', hue='Difference_Flag', hue_order= ['Others(99.51%)', 'DF > 0.05 (0.49%)'], zorder=2, palette=['silver','firebrick'], legend='full')
ax_2.legend(title='')

# Map regression line with scatter plot
plt.plot(x1, predict_y1, color='navy',label='$r^2 = {}$'.format(np.round(r_value_1**2,4))) #from stat
plt.plot(x1, predict_y1, color='navy',label='$y = {} + {}x$'.format(np.round(intercept_1,4), np.round(slope_1,4))) #from stat

plt.text(0,0.9,'$N=12,909$', fontsize=12) #print text on plot
plt.xlabel('TWB1496', fontsize=12)
plt.ylabel('East Asian(EAS, gnomAD)', fontsize=12)
plt.grid(True)

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

#plt.show()
fig_2.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_gnomAD_PASS_DF_r2.png', bbox_inches="tight", dpi=600)

############################################################################### population specific for 'PASS' ########################################################################################

df1 = df[df['Filter']=='PASS(70.68%)']

conditions_eas = [
    (df1['TWB1496'] > 0.01) & (df1['gnomAD_2_Final'] == 0),
    (df1['gnomAD_2_Final'] > 0.01) & (df1['TWB1496'] == 0)]

values_eas = ['TWB', 'AF_eas']

df1['Specific_Flag'] = np.select(conditions_eas, values_eas)
specific = df1['Specific_Flag'].value_counts(normalize=True) * 100
#print(np.round(specific,2))

with open('/work1/viviane1695/tbb1496/updates/1496_EAS_specific.txt', 'w') as f:
    for item in [np.round(specific,2)]:
        f.write("%s\n" % item)

#### Pic 3 ####
select_color = df1.loc[df1['Specific_Flag'] != '0']

x1 = df1['TWB1496']
y1 = df1['gnomAD_2_Final']

slope_1, intercept_1, r_value_1, p_value_1, std_err_1 = stats.linregress(x1,y1)
predict_y1 = intercept_1 + slope_1 * x1

fig_3, ax_3 = plt.subplots()

sns.scatterplot(data=df1, x='TWB1496', y='gnomAD_2_Final', hue='Specific_Flag', hue_order=['0', 'TWB'],palette=['silver', 'lime'], zorder=1, legend=False)

select_color['Specific_Flag'] = select_color['Specific_Flag'].map({'TWB':'TWB specific(0.08%)'})
sns.scatterplot(data=select_color, x='TWB1496', y='gnomAD_2_Final',hue='Specific_Flag', hue_order=['TWB specific(0.08%)'], palette=['lime'], zorder=2, legend='full')
ax_3.legend(title='')

# Map regression line with scatter plot
plt.plot(x1, predict_y1, color='navy',label='$r^2 = {}$'.format(np.round(r_value_1**2,4))) #from stat
plt.plot(x1, predict_y1, color='navy',label='$y = {} + {}x$'.format(np.round(intercept_1,4), np.round(slope_1,4))) #from stat


plt.text(0,0.9,'$N=12,909$', fontsize=12) #print text on plot
plt.xlabel('TWB1496', fontsize=12)
plt.ylabel('East Asian(EAS, gnomAD)', fontsize=12)
plt.grid(True)

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

#plt.show()
fig_3.savefig('/work1/viviane1695/tbb1496/updates/figure/twb1496_gnomAD_PASS_specific.png', bbox_inches="tight", dpi=600)


