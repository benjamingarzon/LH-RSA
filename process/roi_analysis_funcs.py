#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 11:31:00 2020
@author: benjamin.garzon@ki.se
"""
import numpy as np
import os, sys
from nilearn.input_data import NiftiMasker
sys.path.append('/home/xgarzb@GU.GU.SE/Software/')
from PcmPy.matrix import indicator
from sklearn.svm import SVC
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.multiclass import OneVsOneClassifier, OneVsRestClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score, PredefinedSplit
from scipy.optimize import nnls
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score, make_scorer
from glob import glob
import pandas as pd
from sklearn.covariance import ledoit_wolf
from scipy.linalg import svd

minacc = 1.0
MINFEATURES = 500

def project_label(label, hemi, subject, subjects_dir, analysis_dir, labels_dir, overwrite = False):
    commands = []
    commands.append('SUBJECTS_DIR=%s'%(subjects_dir))
    # project label to native surface
    label_orig = os.path.join(labels_dir, '%s.%s.label'%(hemi, label))
    label_dest = os.path.join(analysis_dir, subject, 'label', '%s.%s.label'%(hemi, label))
    label_vol = os.path.join(analysis_dir, subject, 'label', '%s.%s.nii.gz'%(hemi, label))
    
    commands.append('mri_label2label --srclabel %s'%(label_orig) + \
       ' --srcsubject fsaverage' + \
       ' --trgsubject %s'%(subject) + \
       ' --trglabel %s'%(label_dest) + \
       ' --hemi %s --regmethod surface'%(hemi) 
       )
    
    # project to volume
    commands.append( 'mri_label2vol --label %s '%(label_dest) + \
                    ' --subject %s '%(subject) + \
                    ' --hemi %s'%(hemi) + \
                    ' --o %s --regheader T1.mgz --fill-ribbon'%(label_vol)
                    )
    
    command = ';'.join(commands) 
    #print(command)
    if overwrite: 
        print("Projecting label %s to volume"%label)
        os.system(command)
    else:
        if not os.path.exists(label_vol):
            print("Projecting label %s to volume"%label)
            os.system(command)
    return(label_vol)


def second_moment(sequences_sub):
    ncorrect = sequences_sub.shape[0]
    np.sum(sequences_sub.accuracy == 1)
    G_list = []
    G_list.append(np.ones((ncorrect, ncorrect)))
    Xrun = indicator(sequences_sub.run)
    G_list.append(Xrun @ Xrun.T) # Run effect
    Xtrained = 1* (sequences_sub.seq_train.values == 'trained' ).reshape(-1, 1)
    G_list.append(Xtrained @ Xtrained.T) # Trained common
    Xuntrained = 1* (sequences_sub.seq_train.values == 'untrained' ).reshape(-1, 1)
    G_list.append(Xuntrained @ Xuntrained.T) # Untrained common
    
    Xseq1 = indicator(sequences_sub.seq_type*Xtrained.ravel(), positive = True) 
    G_list.append(Xseq1 @ Xseq1.T ) # Difference between two trained sequences
    Xseq2 = indicator(sequences_sub.seq_type*Xuntrained.ravel(), positive = True)
    G_list.append(Xseq2 @ Xseq2.T) # Difference between two untrained sequences
    
    # the last two parameters are the ones of interest
    G_list.append(np.diag(Xtrained.ravel())) # Variability of trained sequences
    G_list.append(np.diag(Xuntrained.ravel())) # Variability of untrained sequences
    
    G = np.dstack(G_list)
    XG = np.hstack([x.reshape(-1, 1) for x in G_list])
    return(XG)

def fit_second_moment(data, XG):
    G_hat = (data @ data.T).ravel()/data.shape[0]
#    return(lsqnonneg(XG, G_hat))
    return(nnls(XG, G_hat))

def second_moment_PCM():
    obs_des = {'cond_vec': cond_vec[i],
               'part_vec': part_vec[i]}
    Y.append(pcm.Dataset(Data[i],obs_descriptors = obs_des))
    N=len(Y)
    G_hat = np.zeros((N,5,5))
    for i in range(N):
        G_hat[i,:,:],_ = pcm.util.est_G_crossval(Y[i].measurements,
                            Y[i].obs_descriptors['cond_vec'],
                            Y[i].obs_descriptors['part_vec'],
                            X=pcm.matrix.indicator(Y[i].obs_descriptors['part_vec']))

    # Build models from the second momement matrices
    M = []
    M.append(pcm.ModelFixed('null',np.eye(5)))
    M.append(pcm.ModelFixed('muscle',modelM[0]))
    M.append(pcm.ModelFixed('natural',modelM[1]))
    M.append(pcm.ModelComponent('muscle+nat',[modelM[0],modelM[1]]))
    M.append(pcm.ModelFree('ceil',5)) # Noise ceiling model 
    T_gr, theta = pcm.inference.fit_model_group(Y, M, fit_scale=True)
    T_cv, theta_cv = pcm.inference.fit_model_group_crossval(Y, M, fit_scale=True)

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
    betas_pre = np.dot(betas[:,:1], np.linalg.inv(cov_ledoit_sqrt))
    return(betas_pre)

def process(label, subject, session, subjects_dir, analysis_dir, labels_dir, 
            sequences, effects_file, residuals_file = None, 
            permutate = False):

    hemi = 'rh' if label.startswith('R_') else 'lh'
    mylabel = project_label(label, hemi, subject, 
                                    subjects_dir, 
                                    analysis_dir, 
                                    labels_dir)
    #print("Opening data for label %s"%label)
    nifti_masker = NiftiMasker(mask_img = mylabel, 
                               memory = os.path.join(analysis_dir, "nilearn_cache"), 
                               memory_level = 1)
    effects = nifti_masker.fit_transform(effects_file)
    if residuals_file:
        # prewhiten the data
        residuals = nifti_masker.fit_transform(residuals_file)
        for run in np.unique(sequences.run):
            effects[sequences.run == run, :] = \
                prewhiten(effects, residuals[sequences.run == run, :])
        
    iscorrect = sequences.accuracy >= minacc
    effects = effects[iscorrect, :]
    targets = sequences.seq_type[iscorrect].array.copy()
    runs = sequences.run[iscorrect]
    if permutate:
        for run in np.unique(runs):
            targets[runs == run] = np.random.permutation(targets[runs == run])
    # second moment
    # print("Computing second moment matrix for label %s"%label)
    XG = second_moment(sequences.loc[iscorrect, :])
    theta, resnorm = fit_second_moment(effects, XG)
    alpha = np.log(theta[7]/theta[6])
    
    # PCM
    # print("Computing SVM classification for label %s"%label)    
    # SVM
    cv_scores = fit_svm(effects, targets, runs)
    return((subject, session, hemi, label, alpha, 
            resnorm, cv_scores, effects.shape[1]))

def gather_results(analysis_dir):
    files = glob(os.path.join(analysis_dir, 'sub-*', 'ses-*', 
                                   'surf', 'roi_scores.csv')
                      )
    results = []
    for file in files:
        results.append(pd.read_csv(file))
        
    results = pd.concat(results)
    return(results)