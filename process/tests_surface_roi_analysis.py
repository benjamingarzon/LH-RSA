#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 13:18:38 2020

@author: xgarzb@GU.GU.SE
"""
import matplotlib.pyplot as plt
from glob import glob
import os
from nilearn.input_data import NiftiMasker
from nilearn.image import resample_to_img
import numpy as np
import pandas as pd
from roi_analysis_funcs import prewhiten, project_label, second_moment, \
    fit_second_moment, second_moment_PCM
from scipy.stats import pearsonr
import seaborn as sns

sns.set_theme(style="white")
cmap = sns.color_palette("coolwarm", as_cmap = True)

WD = '/data/lv0/MotorSkill/'
analysis_dir = '/data/lv0/MotorSkill/fmriprep/analysis/'
subjects_dir = os.path.join(WD, 'fmriprep', 'freesurfer')
labels_dir = '/data/lv0/MotorSkill/labels/fsaverage'

subject = 'sub-lue1204'
session = 1
label = 'R_SPL'
MINACC = 1.0
MINCORRECT = 2
NTYPES = 4

def test_prewhiten():
    
    hemi = 'rh' if label.startswith('R_') else 'lh'
    print("Projecting")
    mylabel = project_label(label, hemi, subject, 
                                    subjects_dir, 
                                    analysis_dir, 
                                    labels_dir, overwrite = False)
    effects_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 'effects.nii.gz')

    resampled_mask = resample_to_img(mylabel, effects_file)

    nifti_masker = NiftiMasker(mask_img = resampled_mask, 
                               memory = os.path.join(analysis_dir, "nilearn_cache"), 
                               memory_level = 1)

    
    sequences_file = os.path.join(WD, 'responses', subject, 
                          'ses-%d'%session, 'sequences.csv')

    run_paths = glob( os.path.join(analysis_dir, subject, 'ses-%d'%session, 
                       'run*'))
    runs = np.array([int(r[-1]) for r in run_paths])
    # find runs
    allsequences = pd.read_csv(sequences_file, sep = ' ')
    sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
    iscorrect = sequences.accuracy >= MINACC
    ncorrect = np.sum(iscorrect)
    effects = nifti_masker.fit_transform(effects_file)
    print(effects.shape)    

    run = 1
    effects_prewhitened = effects.copy()
    residuals_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 
                                     'run%d'%run,
                                     'res4d_LSA_%d.nii.gz'%run)
    residuals = nifti_masker.fit_transform(residuals_file)    
    effects_prewhitened[sequences.run == run, :] = prewhiten(
        effects[sequences.run == run, :], residuals)
    
    plt.scatter(effects_prewhitened[0, :], effects[0, :])

    corr = np.corrcoef(effects)
    plt.figure()
    sns.heatmap(corr, cmap = cmap)
    corr_prew = np.corrcoef(effects_prewhitened)
    plt.figure()
    sns.heatmap(corr_prew, cmap = cmap)
        
    print("Correlation: %f %f"%pearsonr(effects_prewhitened[0, :], 
                                       effects[0, :]))
    
    for run in np.unique(sequences.run):
        print("Opening residuals %d"%run)    
        residuals_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 
                                     'run%d'%run,
                                     'res4d_LSA_%d.nii.gz'%run)
        residuals = nifti_masker.fit_transform(residuals_file)    
        effects[sequences.run == run, :] = prewhiten(
            effects[sequences.run == run, :], residuals)
            

def test_getcorrect():
    sequences_file = os.path.join(WD, 'responses', subject, 
                                  'ses-%d'%session, 'sequences.csv')
    # find runs
    run_paths = glob( os.path.join(analysis_dir, subject, 'ses-%d'%session, 
                   'run*'))
    runs = np.array([int(r[-1]) for r in run_paths])

    allsequences = pd.read_csv(sequences_file, sep = ' ')

    sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
    sequences['iscorrect'] = sequences.accuracy >= MINACC
    # to accept the data of a run, there have to be at least 2 correct instances
    # of each sequence type in it
    valid_types = sequences.groupby(['run', 'seq_type']).iscorrect.sum()
    valid_runs = (valid_types >= MINCORRECT).groupby('run').sum() == NTYPES
    valid_runs = valid_runs.index[valid_runs]
    sequences = sequences.loc[ sequences.run.isin(valid_runs), : ]
    
    ncorrect = np.sum(sequences['iscorrect'])
    return(sequences)

def test_PCM():
    
    hemi = 'rh' if label.startswith('R_') else 'lh'
    print("Projecting")
    mylabel = project_label(label, hemi, subject, 
                                    subjects_dir, 
                                    analysis_dir, 
                                    labels_dir, overwrite = False)
    effects_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 'effects.nii.gz')

    resampled_mask = resample_to_img(mylabel, effects_file)

    nifti_masker = NiftiMasker(mask_img = resampled_mask, 
                               memory = os.path.join(analysis_dir, "nilearn_cache"), 
                               memory_level = 1)

    
    sequences_file = os.path.join(WD, 'responses', subject, 
                          'ses-%d'%session, 'sequences.csv')

    run_paths = glob( os.path.join(analysis_dir, subject, 'ses-%d'%session, 
                       'run*'))
    runs = np.array([int(r[-1]) for r in run_paths])
    # find runs
    allsequences = pd.read_csv(sequences_file, sep = ' ')
    sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
    iscorrect = sequences.accuracy >= MINACC
    ncorrect = np.sum(iscorrect)
    effects = nifti_masker.fit_transform(effects_file)

    # prewhiten the data
    for run in np.unique(sequences.run):
        residuals_file = os.path.join(analysis_dir, subject, 
                                 'ses-%d'%session, 
                                 'run%d'%run,
                                 'res4d_LSA_%d.nii.gz'%run)
        residuals = nifti_masker.fit_transform(residuals_file)    
        effects[sequences.run == run, :], prewhiten_ok = prewhiten(
            effects[sequences.run == run, :], residuals)
    
    iscorrect = sequences.accuracy >= MINACC
    runs = sequences.run[iscorrect]
    
    # PCM
    T_ind, theta, xnobis_trained, xnobis_untrained = second_moment_PCM(sequences.loc[iscorrect, :], effects[iscorrect, :])
    
    return T_ind, theta, xnobis_trained, xnobis_untrained
    # second moment
    # print("Computing second moment matrix for label %s"%label)
    #XG = second_moment(sequences.loc[iscorrect, :])
    #theta, resnorm = fit_second_moment(effects[iscorrect, :], XG)
    #alpha = np.log(theta[7]/theta[6])
    #print(alpha)
#test_prewhiten()
#sequences = test_getcorrect()
T_ind, theta, xnobis_trained, xnobis_untrained = test_PCM()
print(theta[5])
print(T_ind)