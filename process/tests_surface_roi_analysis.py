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
    second_moment_PCM, fit_clf, fit_clf_confusion, compute_distance, \
        get_cov_ledoit, fit_clf_separate
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
from sklearn.decomposition import PCA
from roi_analysis_funcs import svc_ovo, balanced_scorer
from sklearn import manifold
from surface_roi_analysis import do_analysis

sns.set_theme(style="white")
cmap = sns.color_palette("coolwarm", as_cmap = True)

WD = '/data/lv0/MotorSkill/'
analysis_dir = '/data/lv0/MotorSkill/fmriprep/analysis/'
subjects_dir = os.path.join(WD, 'fmriprep', 'freesurfer')
labels_dir = '/data/lv0/MotorSkill/labels/fsaverage'
effects_filename = 'effects.nii.gz'
derivatives_filename = 'derivatives.nii.gz'
subject = 'sub-lue1205'
session = 5
#label = 'R_C1'
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
                                     'ses-%d'%session, effects_filename) 
    derivatives_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, derivatives_filename) 

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
    derivatives = nifti_masker.fit_transform(derivatives_file)

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
                                     'ses-%d'%session, effects_filename)

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
    T_ind, theta, G_hat_trained, G_hat_untrained = second_moment_PCM(sequences.loc[iscorrect, :], effects[iscorrect, :])
    
    return T_ind, theta, G_hat_trained, G_hat_untrained

def test_clf(permutate = False, do_prewhitening = None):
    hemi = 'rh' if label.startswith('R_') else 'lh'
    mylabel = project_label(label, hemi, subject, 
                                    subjects_dir, 
                                    analysis_dir, 
                                    labels_dir, 
                                    overwrite = False)
    
    effects_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, effects_filename)
    derivatives_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, derivatives_filename)

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
    allsequences = allsequences.sort_values(by=['run', 'seq_type'])
    
    sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
    iscorrect = sequences.accuracy >= MINACC
    ncorrect = np.sum(iscorrect)
    effects = nifti_masker.fit_transform(effects_file)
    derivatives = nifti_masker.fit_transform(derivatives_file)

    constant_columns = np.all(effects[1:] == effects[:-1], axis=0)
    effects = effects[:, ~constant_columns]

    if do_prewhitening == 'session':
        # prewhiten the data
        residual_file_list = [os.path.join(analysis_dir, subject, 
                             'ses-%d'%session, 
                             'run%d'%run,
                             'res4d_LSA_%d.nii.gz'%run) \
                      for run in np.unique(sequences.run)]
        cov_ledoit_sqrt = get_cov_ledoit(residual_file_list, nifti_masker, 
                                         constant_columns)
        for run in np.unique(sequences.run):
            effects[sequences.run == run, :], prewhiten_ok = prewhiten(
                effects[sequences.run == run, :], cov_ledoit_sqrt)        

    if do_prewhitening == 'run':
        # prewhiten the data
        for run in np.unique(sequences.run):
            residual_file = os.path.join(analysis_dir, subject, 
                                         'ses-%d'%session, 
                                         'run%d'%run,
                                         'res4d_LSA_%d.nii.gz'%run)
            cov_ledoit_sqrt = get_cov_ledoit([residual_file], nifti_masker, 
                                             constant_columns)
            effects[sequences.run == run, :], prewhiten_ok = prewhiten(
                effects[sequences.run == run, :], cov_ledoit_sqrt)        


    sequences.loc[:, 'iscorrect'] = sequences.accuracy >= MINACC
    # to accept the data of a run, there have to be at least 2 correct instances
    # of each sequence type in it
    valid_types = sequences.groupby(['run', 'seq_type']).iscorrect.sum()
    valid_runs = (valid_types >= MINCORRECT).groupby('run').sum() == NTYPES
    valid_runs = valid_runs.index[valid_runs]
    valid = np.logical_and(sequences.run.isin(valid_runs),
                           sequences.iscorrect)

    #for i in range(4):
    #    fmri_img = index_img(effects_file, np.logical_and(sequences.seq_type == (i + 1), valid))
    #    mean_img(fmri_img).to_filename(os.path.join(analysis_dir, subject, 
    #                                 'ses-%d'%session, 'target-%d.nii.gz'%(i + 1)))

    nvalid = np.sum(valid)
    nn = len(valid)
    nruns = len(valid_runs)
    print("Runs: %d"%nruns)
    # write everything out
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
                
    alldata = pd.concat([pd.Series(nn*[subject]), 
                         pd.Series(nn*[hemi]), 
                         pd.Series(nn*[label]), 
                         pd.Series(nn*[session]), 
                         sequences.run, 
                         sequences.seq_train, 
                         sequences.seq_type, 
                         sequences.iscorrect,
                         valid,
                         pd.DataFrame(effects)], axis = 1)
    alldata.columns = colnames
    sequences = sequences.loc[ valid, : ]
    
    targets = sequences.true_sequence.to_numpy()
    seq_train = sequences.seq_train.to_numpy()
    runs = sequences.run.to_numpy()
    effects = effects[valid, :]
    derivatives = derivatives[valid, :]
    
    targets_perm = np.array(targets, copy = True)

    for run in np.unique(runs):
        targets_perm[runs == run] = \
            np.random.permutation(targets[runs == run])    
            
    if permutate:
        # targets and seq_train are pointing to df
        for run in np.unique(runs):
            myperm = np.random.permutation(range(np.sum(runs == run)))
            targets[runs == run] = targets[runs == run][myperm]
            seq_train[runs == run] = seq_train[runs == run][myperm]
            
    sequences['target'] = sequences.true_sequence   
    df = sequences[['target', 'seq_train']].drop_duplicates()
    #X = np.hstack((effects, derivatives))
    #pca = PCA(n_components = 10).fit(X)
    #Z = pca.transform(X)
    #plt.plot(pca.explained_variance_)
    #plt.show()

    #cv_scores = fit_clf(Z, targets, splits = runs)
    #print(cv_scores)

    #effects = effects[sequences.block == 1, :]
    #targets = targets[sequences.block == 1]
    #runs = runs[sequences.block == 1]
    #effects = np.concatenate((effects, effects))
    #targets = np.concatenate((targets, targets))
    #runs = np.concatenate((runs, runs))
    cosine, cosine_same, cosine_different, cosine_trained_same, \
    cosine_trained_different, cosine_untrained_same, \
    cosine_untrained_different, cosine_trained_untrained, \
    labels = compute_distance(effects, targets, df, splits = runs,
                              dist_type = 'xcosine')

    
    print("cosine: same %f, different %f"%(cosine_same, cosine_different)) 
    xnobis, xnobis_same, xnobis_different, xnobis_trained_same, \
    xnobis_trained_different, xnobis_untrained_same, \
    xnobis_untrained_different, xnobis_trained_untrained, labels = \
        compute_distance(effects, targets, df, splits = runs,  
                         dist_type = 'xnobis')

    print("xnobis: same %f, different %f"%(xnobis_same, xnobis_different))    
        
    cv_scores = fit_clf(effects, targets, splits = runs)
    print(cv_scores)

#    cv_scores = fit_clf(derivatives, targets, splits = runs)
#    print(cv_scores)    
   
    #cv_score, accuracy, accuracy_trained, accuracy_untrained = \
    #    fit_clf_confusion(effects, targets, runs, df)
    clf_acc_trained, clf_acc_untrained,  clf_acc_trained_untrained = \
        fit_clf_separate(effects, targets, runs, df)
    # try grouping
    indices = set([x for x in zip(targets, runs)])
    targets_grouped, runs_grouped = zip(*indices)
    targets_grouped = np.array(targets_grouped)
    runs_grouped = np.array(runs_grouped)
    effects_grouped = np.zeros((len(indices), effects.shape[1]))
    for target, run in  indices:
        effects_grouped[np.logical_and(targets_grouped == target, 
                                       runs_grouped == run), :] = \
            np.mean(effects[np.logical_and(targets == target, runs == run), :], axis = 0)

    #cv_score, accuracy, accuracy_trained, accuracy_untrained = \
    #    fit_clf_confusion(effects_grouped, targets_grouped, runs_grouped, df)
    #cv_scores = fit_clf(effects_grouped, targets_grouped, runs_grouped)        
    
    # by run 
#    cv_scores_list = []
#    for run in np.unique(runs):
#        seq_train[runs == run] = seq_train[runs == run][myperm]
#        X = np.hstack((effects, derivatives))[runs == run, :]

#        pca = PCA(n_components = 4).fit(X)
#        Z = pca.transform(X)

#        cv_scores_list.append(fit_clf(Z,
#                              targets[runs == run], 
#                              splits = sequences.block[runs == run]))
        
#    cv_scores_run = np.mean(np.array(cv_scores_list))
#    print('Scores')
#    print(cv_scores_list, cv_scores_run)  
    
#    stophere
    n_components = 2
    n_neighbors = 10
    methods = {}
    #methods['SE'] = manifold.SpectralEmbedding(n_components=n_components,
    #                                       n_neighbors=n_neighbors)
    methods['t-SNE'] = manifold.TSNE(n_components=n_components, init='pca',
                                 random_state=0)
    #methods['PCA'] = PCA(n_components=n_components)
    markers = {1:'o', 2:'v', 3:'s', 4:'P', 5: '*'}
#    fig = plt.figure(figsize=(15, 8)) 
    # Plot results
    if False: 
        for i, (labelx, method) in enumerate(methods.items()):
            Y = method.fit_transform(effects)
            plt.figure()
            for run in runs:
                plt.scatter(Y[runs== run, 0], Y[runs == run, 1], 
                           c = targets[runs == run], 
                           s = run*20,
                           marker = markers[run],
                           cmap = plt.cm.Spectral)
            plt.title(labelx)
    
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

    if False:    
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
        process_mask[..., :picked_slice] = 0
        process_mask_img = new_img_like(mask_img, process_mask)
    
        fmri_img = index_img(effects_file, valid)

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
    return(cv_scores,            
           ( clf_acc_trained, clf_acc_untrained,  clf_acc_trained_untrained))

def test_do_analysis():
    labels_file='/home/xgarzb@GU.GU.SE/Software/LeftHand/masks/motor_roi_parcels.txt'
    labels_dir='/data/lv0/MotorSkill/labels/fsaverage'
    WD='/data/lv0/MotorSkill/'
    permutate = True
    overwrite_extract = False
    overwrite_scores = True
    output_data = False 
    do_prewhitening = True
    suffix = 'roi-mask-perm'
    effects_name = 'effects.nii.gz' 
    num_cores = 1
    
    do_analysis(WD, permutate, overwrite_extract, 
                overwrite_scores, output_data, 
                do_prewhitening, suffix, labels_file, 
                labels_dir, effects_name, num_cores)
#test_prewhiten()
#sequences = test_getcorrect()
#T_ind, theta, G_hat_trained, G_hat_untrained = test_PCM()

#test_PCM()
#for elem in test_clf(permutate = False, do_prewhitening = 'session'):
#    print(elem)
#for elem in test_clf(permutate = True, do_prewhitening = 'session'):
#    print(elem)

test_do_analysis()
