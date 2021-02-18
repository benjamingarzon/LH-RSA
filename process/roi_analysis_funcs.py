#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 11:31:00 2020
@author: benjamin.garzon@ki.se
"""
import matplotlib.pyplot as plt
import time
import numpy as np
import os, sys
import pandas as pd
from glob import glob
from nilearn.input_data import NiftiMasker
from nilearn.image import resample_to_img
from sklearn.svm import SVC
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.multiclass import OneVsOneClassifier, OneVsRestClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score, PredefinedSplit
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score, make_scorer
from sklearn.covariance import ledoit_wolf
from scipy.optimize import nnls
from scipy.linalg import svd
sys.path.append('/home/xgarzb@GU.GU.SE/Software/')
import PcmPy as pcm

MINACC = 1.0
MINFEATURES = 'all'
MINCORRECT = 3
NTYPES = 4

def project_label(label, hemi, subject, subjects_dir, analysis_dir, labels_dir, 
                  overwrite = False):
    
    ribbon_mask = os.path.join(analysis_dir, subject, 'label', 
                               'ribbon_mask.nii.gz')

    commands = []
    commands.append('SUBJECTS_DIR=%s'%(subjects_dir))
    # project label to native surface
    label_orig = os.path.join(labels_dir, '%s.%s.label'%(hemi, label))
    label_dest = os.path.join(analysis_dir, subject, 'label', '%s.%s.label'%(hemi, label))
    label_vol = os.path.join(analysis_dir, subject, 'label', '%s.%s.nii.gz'%(hemi, label))
    struct_img = os.path.join(subjects_dir, subject, 'mri', 'T1.mgz')
    
    commands.append('mri_label2label --srclabel %s'%(label_orig) + \
       ' --srcsubject fsaverage' + \
       ' --trgsubject %s'%(subject) + \
       ' --trglabel %s'%(label_dest) + \
       ' --hemi %s --regmethod surface'%(hemi) 
       )
    
    # project to volume
    commands.append( 'mri_label2vol --label %s '%(label_dest) + \
                    ' --subject %s '%(subject) + \
                    ' --hemi %s '%(hemi) + \
                    ' --o %s --regheader %s '%(label_vol, struct_img) + \
                    '--temp %s ' %(struct_img) + \
                    '--proj frac 0.1 0.9 0.1' 
                    ) #fill-ribbon
        
    commands.append('fslmaths %s -dilF -ero -mas %s %s -odt char'%(label_vol, 
                                                                   ribbon_mask, 
                                                                   label_vol))
    command = ';'.join(commands) 
    #print(command)
    if overwrite: 
        print("Projecting label %s to volume"%label)
        os.system(command)
    else:
        if not os.path.exists(label_vol):
            print("Projecting label %s to volume"%label)
            try:
                if os.system(command) != 0:
                    print(command)
                    raise Exception(command)
            except:
                print("Command failed")
    return(label_vol)


def second_moment(sequences_sub):
    ncorrect = sequences_sub.shape[0]
    G_list = []
    G_list.append(np.ones((ncorrect, ncorrect)))
    Xrun = pcm.matrix.indicator(sequences_sub.run)
    G_list.append(Xrun @ Xrun.T) # Run effect
    Xtrained = 1*(sequences_sub.seq_train.values == 'trained').reshape(-1, 1)
    G_list.append(Xtrained @ Xtrained.T) # Trained common
    Xuntrained = 1*(sequences_sub.seq_train.values == 'untrained').reshape(-1, 1)
    G_list.append(Xuntrained @ Xuntrained.T) # Untrained common
    
    Xseq1 = pcm.matrix.indicator(sequences_sub.seq_type*Xtrained.ravel(), 
                                 positive = True) 
    G_list.append(Xseq1 @ Xseq1.T ) # Difference between two trained sequences
    Xseq2 = pcm.matrix.indicator(sequences_sub.seq_type*Xuntrained.ravel(), 
                                 positive = True)
    G_list.append(Xseq2 @ Xseq2.T) # Difference between two untrained sequences
    
    # the last two parameters are the ones of interest
    G_list.append(np.diag(Xtrained.ravel())) # Variability of trained sequences
    G_list.append(np.diag(Xuntrained.ravel())) # Variability of untrained sequences
    
    G = np.dstack(G_list)
    XG = np.hstack([x.reshape(-1, 1) for x in G_list])
    return(XG)

def fit_second_moment(data, XG):
    G_hat = (data @ data.T).ravel()/data.shape[0]
    try:
        theta, resnorm = nnls(XG, G_hat)
    except RuntimeError:
        theta = np.zeros(XG.shape[1])
        resnorm = -1
    return(theta, resnorm)

def second_moment_PCM(sequences_sub, data):
    
    Ntypes = len(np.unique(sequences_sub.seq_type))
    runs = np.unique(sequences_sub.run)
    Nruns = len(runs)
    obs_des = {'cond_vec': sequences_sub.seq_type,
               'part_vec': sequences_sub.run}
    Y = pcm.Dataset(data, obs_descriptors = obs_des)
    
    G_hat,_ = pcm.util.est_G_crossval(Y.measurements,
                            Y.obs_descriptors['cond_vec'],
                            Y.obs_descriptors['part_vec'],
                            X = pcm.matrix.indicator(
                                Y.obs_descriptors['part_vec']))
    
    vmin = G_hat.min()
    vmax = G_hat.max()
    #plt.imshow(G_hat,vmin=vmin,vmax=vmax)

    # define models 
    df = sequences_sub[['seq_type', 'seq_train']].drop_duplicates()
    mapping = pd.Series(df.seq_train.values, index=df.seq_type) #.to_dict()
    xnobis = np.diag(G_hat)
    xnobis_trained = np.mean(xnobis[mapping == 'trained'])
    xnobis_untrained = np.mean(xnobis[mapping == 'untrained'])
    
    G_list = []    
    G_list.append(np.ones((Ntypes, Ntypes)))
    Xtrained = 1* (mapping == 'trained').values.reshape(-1, 1)
    G_list.append(Xtrained @ Xtrained.T) # Trained common
    Xuntrained = 1* (mapping == 'untrained' ).values.reshape(-1, 1)
    G_list.append(Xuntrained @ Xuntrained.T) # Untrained common
    
    trained_diag = np.diag(Xtrained.ravel())
    untrained_diag = np.diag(Xuntrained.ravel())

    Xseq1 = G_list[1] - trained_diag  
    G_list.append(Xseq1) # Difference between two trained sequences
    Xseq2 = G_list[2] - untrained_diag
    G_list.append(Xseq2) # Difference between two untrained sequences

    # the last two parameters are the ones of interest
    G_list.append(trained_diag) # Variability of trained sequences
    G_list.append(untrained_diag) # Variability of untrained sequences

    M = pcm.ModelComponent('model', G_list)
    
    #for i in range(len(G_list)):
    #    print(i)
    #    plt.subplot(3,3,i+1)
    #    plt.imshow(G_list[i])
   
    T_ind, theta = \
        pcm.inference.fit_model_individ(Y, M, fit_scale=True)
    # last two theta are scaling and noise
    return(T_ind, theta[0].ravel(), xnobis_trained, xnobis_untrained)


svc_ovo = OneVsOneClassifier(Pipeline([
        ('standardscaler', StandardScaler()),
        ('anova', SelectKBest(f_classif, k = MINFEATURES)),
        ('svc', SVC(kernel='linear'))
    ]))
    
svc_ova = OneVsRestClassifier(Pipeline([
    ('standardscaler', StandardScaler()),
    ('anova', SelectKBest(f_classif, k = MINFEATURES)),
    ('svc', SVC(kernel='linear'))
]))

balanced_scorer = make_scorer(balanced_accuracy_score)

def fit_svm(X, y, runs):
    cv_scores = cross_val_score(svc_ovo, X, y, cv = PredefinedSplit(runs), 
                                verbose = 0, 
                                scoring = balanced_scorer)
    return(cv_scores.mean())
    #print('OvO:', cv_scores_ovo.mean())
    
def prewhiten(betas, residuals):
    cov_ledoit, _ = ledoit_wolf(residuals)    
    Uc, Dc, Vhc = svd(cov_ledoit, full_matrices = False)
    cov_ledoit_sqrt = np.dot(Uc * np.sqrt(Dc), Vhc)
    prewhiten_ok = True
    try:
        betas_prewhitened = np.dot(betas, np.linalg.inv(cov_ledoit_sqrt))
    except np.linalg.LinAlgError:
        betas_prewhitened = betas #if there are problems use original
        prewhiten_ok = False
    return(betas_prewhitened, prewhiten_ok)

def process(label, subject, session, subjects_dir, analysis_dir, labels_dir, 
            sequences, effects_file, do_prewhitening = False, 
            permutate = False):
            
    hemi = 'rh' if label.startswith('R_') else 'lh'
    mylabel = project_label(label, hemi, subject, 
                                    subjects_dir, 
                                    analysis_dir, 
                                    labels_dir)

    #print("Opening data for label %s"%mylabel)
    resampled_mask = resample_to_img(mylabel, effects_file, 
                                     interpolation = 'nearest')
    nifti_masker = NiftiMasker(mask_img = resampled_mask, 
                               memory = os.path.join(analysis_dir, "nilearn_cache"), 
                               memory_level = 1)
    try:
        effects = nifti_masker.fit_transform(effects_file)
    except ValueError:
        print("Could no open data for %s. Input data: %s"%(mylabel, effects_file))
   
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
    
    sequences = sequences.loc[ valid, : ]
    
    targets = sequences.seq_type.array.copy()
    runs = sequences.run
    effects = effects[valid, :]
    
    if permutate:
        for run in np.unique(runs):
            targets[runs == run] = np.random.permutation(targets[runs == run])
    # second moment
    # print("Computing second moment matrix for label %s"%label)
    XG = second_moment(sequences)
    theta, resnorm = fit_second_moment(effects, XG)
    T_ind, theta_PCM, xnobis_trained, xnobis_untrained = \
        second_moment_PCM(sequences, effects)

    alpha_trained = np.log(theta[6])
    alpha_untrained = np.log(theta[7])

    alpha_trained_PCM = np.log(theta_PCM[5])
    alpha_untrained_PCM = np.log(theta_PCM[6])
    
    # PCM
    # print("Computing SVM classification for label %s"%label)    
    # SVM
    cv_scores = fit_svm(effects[sequences.iscorrect, :], targets, runs)
    return((subject, session, hemi, label, resnorm,
            alpha_trained, 
            alpha_untrained, 
            alpha_trained_PCM, 
            alpha_untrained_PCM, 
            xnobis_trained,
            xnobis_untrained,
            cv_scores, effects.shape[1]))

def gather_results(analysis_dir):
    files = glob(os.path.join(analysis_dir, 'sub-*', 'ses-*', 
                                   'surf', 'roi_scores.csv')
                      )
    results = []
    for file in files:
        results.append(pd.read_csv(file))
        
    results = pd.concat(results)
    return(results)