#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 11:17:22 2021

@author: xgarzb@GU.GU.SE
"""
import pandas as pd
from scipy.stats import iqr
import seaborn as sns
import os
import matplotlib.pyplot as plt

WD = '/data/lv0/MotorSkill/'

def spot_outliers_IQR(x):
    Q1=x.quantile(0.25)
    Q3=x.quantile(0.75)
    IQR=Q3-Q1
    thr = (Q1-1.5*IQR)
    return x < thr, thr

results_df = pd.read_csv(os.path.join(WD, 'fmriprep/fmriprep/reg_stats.txt'))


sns.boxplot(data = results_df, x = 'isrejected', y = 'MI')
plt.figure()
sns.boxplot(data = results_df, x = 'wave', y = 'MI')
plt.figure()
sns.histplot(data = results_df, x = 'MI')


outliers, thr = spot_outliers_IQR(results_df.MI)
results_df.sort_values('MI').head(20)
maxrows = 20
i = 0

for index, row in results_df.sort_values('MI').iterrows():
    print(row)
    os.system(
        'firefox file:///data/lv0/MotorSkill/fmriprep/fmriprep/{}/figures/{}_ses-{}_task-sequence_run-{}_desc-coreg_bold.svg'.format(
            row['subject'], row['subject'], row['session'], row['run']))
    print(row['subject'], row['session'], row['run'])
    if i == maxrows:
        break
    i = i + 1