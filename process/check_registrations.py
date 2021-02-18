#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# values above 0.15 seem ok
"""
Created on Thu Nov 12 11:31:00 2020
@author: benjamin.garzon@ki.se
"""

import multiprocessing
from joblib import Parallel, delayed
from scipy.stats import pearsonr
from nilearn.image import resample_to_img, load_img, mean_img
from nilearn.input_data import NiftiMasker
from nilearn.masking import intersect_masks
from sklearn.metrics import normalized_mutual_info_score, mutual_info_score
from sklearn.feature_selection import mutual_info_regression
import pandas as pd
import os
import numpy as np
import pandas as pd
from glob import glob
import seaborn as sns
import matplotlib.pyplot as plt
from nilearn.plotting import plot_anat, plot_epi
    
# project the labels to volume
#if __name__ == "__main__":
    
WD = '/data/lv0/MotorSkill/'
analysis_dir = '/data/lv0/MotorSkill/fmriprep/fmriprep/'
paths = glob(os.path.join(analysis_dir, 'sub-*', 'ses-*', 'func',
                              '*T1w_desc-preproc_bold.nii.gz'))
    
# iterate across paths
results = []
for path in paths:
    pathels = os.path.normpath(path).split(os.sep)
    subject = pathels[-4]
    session = int(pathels[-3].split('-')[1])
    run = int(pathels[-1].split('_')[3].split('-')[1])
    
    struct_file = os.path.join(analysis_dir, subject, 'anat',
                             '%s_desc-preproc_T1w.nii.gz'%subject)
    struct_mask_file = os.path.join(analysis_dir, subject, 'anat',
                     '%s_desc-brain_mask.nii.gz'%subject)

    bold_mask_file = glob(os.path.join(analysis_dir, subject, 
                                  'ses-%d'%session, 
                                  'func',
                                  '%s_ses-%d_task-sequence_run-%d_space-T1w_desc-brain_mask.nii.gz'%(subject, session, run)))        
    bold_file = path            
    struct_img = resample_to_img(struct_file, bold_file)
    struct_mask_img = resample_to_img(struct_mask_file, bold_mask_file)
    bold_img = mean_img(bold_file)
    bold_mask_img = load_img(bold_mask_file)
    mask = intersect_masks((struct_mask_img, bold_mask_img))
    masker = NiftiMasker(mask)
    x = masker.fit_transform(struct_img).ravel()
    y = masker.fit_transform(bold_img).ravel()
    ind = np.logical_and(x>0, y>0)
    x = x[ind]
    y = y[ind]
    MI = mutual_info_regression(x.reshape(-1, 1), y)[0]
    cor = pearsonr(x, y)[0]
    res = [subject, session, run, MI, cor]
    print(res)
    results.append(res)
    #plot_epi(struct_img, cut_coords = (0, 0, 0))
    #plot_epi(bold_img, cut_coords = (0, 0, 0))
    #plot_epi(mask, cut_coords = (0, 0, 0))
    #sns.scatterplot(np.log(x), np.log(y))
    #plt.clf()
results_df = pd.DataFrame(results, columns = ('subject','session','run', 'MI', 'cor'))
sns.scatterplot(data = results_df, x = 'MI', y = 'cor')
plt.hist(results_df.MI, 100)

# compare with rejected bbregister from fmriprep
rejected = pd.read_csv(os.path.join(WD, 
                                 'fmriprep/fmriprep/rejected_regs.txt'), 
                                 header = None, 
                                 names = ['subject', 'session', 'run'],
                                 sep = ' ').apply(lambda x: x.str.split('-').str[1].astype(int) if x.name in ['session', 'run'] else x)
rejected['isrejected'] = 1
results_df = pd.merge(results_df, 
                      rejected, 
                      how = "left", 
                      on = ['subject', 'session', 'run'])
results_df['isrejected'] = results_df['isrejected'].isnull().apply(np.invert)

print(results_df.head())
sns.boxplot(data = results_df, x = 'isrejected', y = 'MI')
