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
from roi_analysis_funcs import prewhiten, project_label, \
    second_moment_PCM, fit_svm, fit_svm_confusion
from scipy.stats import pearsonr
import seaborn as sns
from nilearn.image import index_img
import nilearn.decoding
from nilearn.image.image import mean_img, math_img
from nilearn.plotting import plot_epi, show
from nilearn.masking import compute_epi_mask
from nilearn.plotting import plot_roi
from nilearn.image import new_img_like, load_img, get_data
from sklearn.model_selection import cross_val_score, PredefinedSplit, \
GridSearchCV
from roi_analysis_funcs import svc_ovo, balanced_scorer

sns.set_theme(style="white")
cmap = sns.color_palette("coolwarm", as_cmap = True)

WD = '/data/lv0/MotorSkill/'
analysis_dir = '/data/lv0/MotorSkill/fmriprep/analysis/'
subjects_dir = os.path.join(WD, 'fmriprep', 'freesurfer')
labels_dir = '/data/lv0/MotorSkill/labels/fsaverage'

subject = 'sub-lue1101'
session = 2
label = 'R_C1'
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

def test_svm(permutate = True, do_prewhitening = True):
    hemi = 'rh' if label.startswith('R_') else 'lh'
    mylabel = project_label(label, hemi, subject, 
                                    subjects_dir, 
                                    analysis_dir, 
                                    labels_dir, 
                                    overwrite = False)
    
    effects_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 'effects.nii.gz')

    searchlight_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 'searchlight.nii.gz')
    
    resampled_mask = resample_to_img(mylabel, effects_file, 
                                     interpolation = 'nearest')

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
    
    # reorder in time to sync with data!
    allsequences = allsequences.sort_values(by=['run', 'onset'])
    
    sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
    iscorrect = sequences.accuracy >= MINACC
    ncorrect = np.sum(iscorrect)
    effects = nifti_masker.fit_transform(effects_file)

    if do_prewhitening:
        # prewhiten the data
        for run in np.unique(sequences.run):
            residuals_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 
                                     'run%d'%run,
                                     'res4d_LSA_%d.nii.gz'%run)
            residuals = nifti_masker.fit_transform(residuals_file)    
            effects[sequences.run == run, :], prewhiten_ok = prewhiten(
                effects[sequences.run == run, :], residuals)

    sequences.loc[:, 'iscorrect'] = sequences.accuracy >= MINACC
    # to accept the data of a run, there have to be at least 2 correct instances
    # of each sequence type in it
    valid_types = sequences.groupby(['run', 'seq_type']).iscorrect.sum()
    valid_runs = (valid_types >= MINCORRECT).groupby('run').sum() == NTYPES
    valid_runs = valid_runs.index[valid_runs]
    valid = np.logical_and(sequences.run.isin(valid_runs),
                           sequences.iscorrect)

    for i in range(4):
        fmri_img = index_img(effects_file, np.logical_and(sequences.seq_type == (i + 1), valid))
        mean_img(fmri_img).to_filename(os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 'target-%d.nii.gz'%(i + 1)))

    fmri_img = index_img(effects_file, np.logical_and(sequences.seq_type == (i + 1), valid))

    # write everything out
    nn = len(valid)
    colnames = ['subject',
                'hemi',
                'label',
                'session',
                'run',
                'seq_train',
                'seq_type',
                'iscorrect', 
                'valid'] + \
        ['feature_%0.4d'%(x) for x in range(effects.shape[1])]
                
    alldata = pd.concat([pd.Series(nn*subject), 
                         pd.Series(nn*hemi), 
                         pd.Series(nn*label), 
                         pd.Series(nn*session), 
                         sequences.run, 
                         sequences.seq_train, 
                         sequences.seq_type, 
                         sequences.iscorrect,
                         valid,
                         pd.DataFrame(effects)], axis = 1)
    alldata.columns = colnames
    sequences = sequences.loc[ valid, : ]
    
    targets = sequences.seq_type.array.copy()
    runs = sequences.run
    effects = effects[valid, :]
    
    if permutate:
        for run in np.unique(runs):
            targets[runs == run] = np.random.permutation(targets[runs == run])
    df = sequences[['seq_type', 'seq_train']].drop_duplicates()
    cv_score, accuracy, accuracy_trained_accuracy_untrained = \
        fit_svm_confusion(effects, targets, runs, df)
    cv_scores = fit_svm(effects, targets, runs)
    
    from sklearn import manifold
    n_components = 2
    n_neighbors = 10
    methods = {}
    #methods['SE'] = manifold.SpectralEmbedding(n_components=n_components,
    #                                       n_neighbors=n_neighbors)
    methods['t-SNE'] = manifold.TSNE(n_components=n_components, init='pca',
                                 random_state=0)

#    fig = plt.figure(figsize=(15, 8)) 
    # Plot results
    for i, (labelx, method) in enumerate(methods.items()):
        Y = method.fit_transform(effects)
        plt.scatter(Y[:, 0], Y[:, 1], 
                   c = targets,
                   s = runs*20, 
                   cmap = plt.cm.Spectral)

    #print(runs.value_counts())
#    cv_scores = cross_val_score(svc_ovo, X, targets, cv = PredefinedSplit(runs), 
#                                verbose = 0, 
#                                scoring = balanced_scorer)
    mean_signal_trained = \
        np.mean(
            np.mean(effects[sequences.seq_train.to_numpy() == 'trained' , :], 
                    axis = 0))
    mean_signal_untrained = \
        np.mean(
            np.mean(effects[sequences.seq_train.to_numpy() == 'untrained' , :], 
            axis = 0))
    print(mean_signal_trained, mean_signal_untrained)

    if True:    
        # do searchlight analysis
            
        # Compute the mean EPI: we do the mean along the axis 3, which is time
        mask_img = math_img("np.std(img1, axis = -1)> 10", img1 = effects_file)
        
        # Visualize it as an ROI
        mean_effects = mean_img(effects_file)
        plot_roi(mask_img, mean_effects)
        show()
    
        process_mask = get_data(mask_img).astype(np.int)
        picked_slice = 52
    #    process_mask[..., (picked_slice + 3):] = 0
        #process_mask[..., :picked_slice] = 0
        process_mask_img = new_img_like(mask_img, process_mask)
    
        fmri_img = index_img(effects_file, valid)
    #    cv_scores = cross_val_score(cl, X, y, cv = PredefinedSplit(runs), 
    #                            verbose = 0, 
    #                            scoring = balanced_scorer)
        # y, session = y[condition_mask], session[condition_mask]
        # # The radius is the one of the Searchlight sphere that will scan the volume
        searchlight = nilearn.decoding.SearchLight(
        mask_img,
        estimator = svc_ovo,
        scoring = balanced_scorer,
        process_mask_img=process_mask_img,
        radius=6, n_jobs=10,
        verbose=1, cv=PredefinedSplit(runs))
    
        searchlight.fit(fmri_img, targets)
        searchlight_img = new_img_like(mask_img, searchlight.scores_)
        searchlight_img.to_filename(searchlight_file)
    return(cv_score, accuracy, accuracy_trained_accuracy_untrained)
#test_prewhiten()
#sequences = test_getcorrect()
#T_ind, theta, xnobis_trained, xnobis_untrained = test_PCM()


print(test_svm(permutate = False, do_prewhitening = False))