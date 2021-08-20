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
from sklearn.feature_selection import SelectKBest, f_classif, VarianceThreshold
from sklearn.multiclass import OneVsOneClassifier, OneVsRestClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score, LeaveOneGroupOut, \
    GridSearchCV, KFold
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score, make_scorer
from sklearn.covariance import ledoit_wolf
from scipy.optimize import nnls
from scipy.linalg import svd
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_validate, cross_val_predict, LeaveOneOut
from scipy.spatial.distance import pdist, cdist, cosine, correlation
from itertools import combinations, product
import random

sys.path.append('/home/xgarzb@GU.GU.SE/Software/')
import PcmPy as pcm

MINACC = 1.0
MINFEATURES = 'all'
MINCORRECT = 3
NTYPES = 4

def project_label(label, hemi, subject, subjects_dir, analysis_dir, labels_dir, 
                  overwrite = True):
    
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
                            X = pcm.matrix.indicator(Y.obs_descriptors['part_vec']))
    
    vmin = G_hat.min()
    vmax = G_hat.max()
    #plt.imshow(G_hat,vmin=vmin,vmax=vmax)

    # define models 
    df = sequences_sub[['seq_type', 'seq_train']].drop_duplicates()
    mapping = pd.Series(df.seq_train.values, index=df.seq_type) #.to_dict()
    G_hat_d = np.diag(G_hat)
    G_hat_trained = np.mean(G_hat_d[mapping == 'trained'])
    G_hat_untrained = np.mean(G_hat_d[mapping == 'untrained'])
    
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
    return(T_ind, theta[0].ravel(), G_hat_trained, G_hat_untrained)


param_grid_simple = {
    "svc__C": [np.exp(i) for i in np.arange(-12, -1, 1)]
    }  

param_grid0 = {
    "estimator__svc__C": [np.exp(i) for i in np.arange(-12, -1, 1)]
    }  

param_grid1 = {
    "estimator__svc__C": [np.exp(i) for i in np.arange(-12, -1, 0.5)],
    "estimator__svc__gamma": [np.exp(i) for i in np.arange(-3, 1, 0.5)]
    }  

param_grid_pca = {
    "estimator__svc__C": [np.exp(i) for i in np.arange(-12, -1, 0.5)]
  #  "estimator__pca__n_components": [10, 20, 30, 40, 50, 60]
    }  

param_grid_pls = {
    "estimator__pls__n_components": [10, 20, 30, 40, 50]
    }


inner_cv = KFold(n_splits=5, shuffle=True)

outer_cv = KFold(n_splits=5, shuffle=True)
  
balanced_scorer = make_scorer(balanced_accuracy_score)

svc0 = Pipeline([
        ('variancethreshold', VarianceThreshold()),
#        ('pca', PCA(n_components = 5)),
        ('standardscaler', StandardScaler()),
        ('anova', SelectKBest(f_classif, k = MINFEATURES)), # 
        ('svc', SVC(kernel='linear'))
    ])

svc1 = Pipeline([
        ('variancethreshold', VarianceThreshold()),
        ('standardscaler', StandardScaler()),
        ('anova', SelectKBest(f_classif, k = MINFEATURES)),
        ('svc', SVC(kernel='rbf'))
    ])

pls = Pipeline([
        ('variancethreshold', VarianceThreshold()),
        ('standardscaler', StandardScaler()),
        ('anova', SelectKBest(f_classif, k = MINFEATURES)),
        ('pls', PLSRegression(n_components = 10))
    ])


svc_ova = OneVsRestClassifier(svc0)

svc_ovo = OneVsOneClassifier(svc0)

svc_ova_nl = OneVsRestClassifier(svc1)

svc_ovo_nl = OneVsOneClassifier(svc1)

pls_ovo = OneVsOneClassifier(pls)

balanced_scorer = make_scorer(balanced_accuracy_score)

def partial_dist(X, y, splits, label1, label2, splitp, distance = cosine):
    
    if label1 == label2:
        indicesp = np.where(np.logical_and(y == label1, splits == splitp))[0]
        index_pairsp = list(combinations(indicesp, 2))

    else:
        indices1p = np.where(np.logical_and(y == label1, splits == splitp))[0]
        indices2p = np.where(np.logical_and(y == label2, splits == splitp))[0]
        index_pairsp = list(product(indices1p, indices2p))

    dists = [distance(X[index1, :], X[index2, :]) \
                 for index1, index2 in index_pairsp ]
    print(dists)
    return np.mean(dists)

def partial_xnobis(X, y, splits, label1, label2, splitp, splitn, n_sample = None):

    if label1 == label2:
        indicesp = np.where(np.logical_and(y == label1, splits == splitp))[0]
        indicesn = np.where(np.logical_and(y == label1, splits == splitn))[0]
        index_pairsp = list(combinations(indicesp, 2))
        index_pairsn = list(combinations(indicesn, 2))

    else:
        indices1p = np.where(np.logical_and(y == label1, splits == splitp))[0]
        indices1n = np.where(np.logical_and(y == label1, splits == splitn))[0]
        indices2p = np.where(np.logical_and(y == label2, splits == splitp))[0]
        indices2n = np.where(np.logical_and(y == label2, splits == splitn))[0]
        index_pairsp = list(product(indices1p, indices2p))
        index_pairsn = list(product(indices1n, indices2n))

    index_prod = list(product(index_pairsp, index_pairsn))
    if n_sample is not None:
        if n_sample < len(index_prod):
            index_prod = random.sample(index_prod, k = n_sample)
                
    dists = [np.dot(X[indexp[0], :] -  X[indexp[1], :], 
                    X[indexn[0], :] -  X[indexn[1], :]) \
             for indexp, indexn in index_prod ]
           
    return np.mean(dists)/X.shape[1]

  
def compute_distance(X, y, df, splits, dist_type = 'xnobis'):
    
    labels = np.unique(y)
    n_labels = len(labels)
    usplits = np.unique(splits)

    if dist_type == 'xnobis':

        dist = np.zeros((n_labels, n_labels, len(usplits), len(usplits)-1))
    
        for j, splitp in enumerate(usplits):
            for k, splitn in enumerate(usplits[usplits != splitp]):
                for i1, label1 in enumerate(labels):
                    for i2, label2 in enumerate(labels):
                            dist[i1, i2, j, k] = partial_xnobis(X, y, splits, 
                                                                  label1, label2, 
                                                                  splitp, splitn)
        dist = np.mean(dist, axis = 3) # across training splits 
 
    else:

        dist = np.zeros((n_labels, n_labels, len(usplits)))
        distance = correlation if dist_type == 'correlation' else cosine
        for j, splitp in enumerate(usplits):
            for i1, label1 in enumerate(labels):
                for i2, label2 in enumerate(labels):
                            dist[i1, i2, j] = partial_dist(X, y, 
                                                             splits, 
                                                             label1, label2, 
                                                             splitp, 
                                                             distance)
                        
    dist = np.mean(dist, axis = 2) # across validation splits
    
    plt.figure()
    plt.imshow(dist, cmap = 'jet')
    dist_same = np.mean(np.diag(dist))
    dist_different = np.mean(dist[np.eye(dist.shape[0])==0])
    
    mapping = pd.Series(df.seq_train.values, index=df.seq_type) 
 
    dist_trained = dist[mapping == 'trained', :][:, mapping == 'trained']
    dist_untrained = dist[mapping == 'untrained', :][:, mapping == 'untrained']
    dist_trained_untrained = np.vstack(
        (dist[mapping == 'trained', :][:, mapping == 'untrained'],
         dist[mapping == 'untrained', :][:, mapping == 'trained'])
        )

    elem_trained = np.eye(dist_trained.shape[0]) == 0
    elem_untrained = np.eye(dist_untrained.shape[0]) == 0
    
    dist_trained_same = np.mean(np.diag(dist_trained))
    dist_trained_different = np.mean(dist_trained[elem_trained])
    dist_untrained_same = np.mean(np.diag(dist_untrained))
    dist_untrained_different = np.mean(dist_untrained[elem_untrained])
    dist_trained_untrained = np.mean(dist_trained_untrained)
         
    return(dist, 
           dist_same, 
           dist_different, 
           dist_trained_same, 
           dist_trained_different, 
           dist_untrained_same,
           dist_untrained_different, 
           dist_trained_untrained,
           labels)

def fit_clf(X, y, splits = None, cl = svc_ovo):
    
    if True:
        # with grid_search
        cl_cv = GridSearchCV(estimator = cl, 
                              param_grid = param_grid0, 
                              cv = LeaveOneGroupOut(), 
                              verbose = 0, 
                              scoring = balanced_scorer)

        cv_scores = cross_val_score(cl_cv, X, y, cv = LeaveOneGroupOut(), 
                                    verbose = 0, groups = splits, 
                                    scoring = balanced_scorer,
                                    fit_params={"groups": splits})

    else:
        cv_scores = cross_val_score(cl, X, y, cv = LeaveOneGroupOut(), 
                                    verbose = 0,  groups = splits, 
                                    scoring = balanced_scorer)

    if True:
        cl_cv.fit(X, y, groups = splits)
        print(cl_cv.best_params_)
        #print(cv_scores)
        plt.figure()
        plt.plot(param_grid0['estimator__svc__C'], 
                 cl_cv.cv_results_['mean_test_score'])
    # plt.show()
    return(cv_scores.mean())


def fit_clf_confusion(X, y, splits, df, cl = svc_ovo):
    mapping = pd.Series(df.seq_train.values, index=df.seq_type) 
    
    # cl_cv = GridSearchCV(estimator = cl, 
    #                   param_grid=param_grid, 
    #                   cv=LeaveOneGroupOut(runs), verbose = 0, 
    #                   scoring = balanced_scorer)
    y_pred = cross_val_predict(cl, X, y, cv = LeaveOneGroupOut(), 
                               groups = splits, verbose = 0)
    
    conf_mat = confusion_matrix(y, y_pred)
    accuracies = conf_mat.diagonal()/np.sum(conf_mat, axis = 1)
    accuracy_trained = np.mean(accuracies[mapping == 'trained'])
    accuracy_untrained = np.mean(accuracies[mapping == 'untrained'])
    
#    cv_results = cross_validate(cl, X, y, cv = LeaveOneGroupOut(runs),
#                                verbose = 1,
#                                scoring = confusion_matrix_scorer)
#    cv_scores = cross_val_score(cl, X, y, cv = LeaveOneGroupOut(runs), 
#                                verbose = 0, 
#                                scoring = balanced_scorer)

    return(accuracies, accuracies.mean(), accuracy_trained, accuracy_untrained)   

def fit_clf_separate(X, y, splits, df, cl = svc0):
    
    cl_cv = GridSearchCV(estimator = cl, 
                         param_grid=param_grid_simple, 
                         cv=LeaveOneGroupOut(),
                         verbose = 0, 
                         scoring = balanced_scorer)
    
    mapping = pd.Series(df.seq_train.values, index=df.seq_type) 
 
    labels = np.unique(y)
    
    trained_scores = []
    untrained_scores = []
    trained_untrained_scores = []
    
    for labelx in labels:
        for labely in labels:
            if labelx != labely:
                indices = np.logical_or(y == labelx, y == labely)
                cv_scores = cross_val_score(cl_cv, X[indices, :], y[indices], 
                                        cv = LeaveOneGroupOut(), 
                                        verbose = 0,  groups = splits[indices], 
                                        scoring = balanced_scorer,
                                        fit_params={"groups": splits[indices]})
                if mapping[labelx] != mapping[labely]:
                    trained_untrained_scores.extend(cv_scores.tolist())                    
                elif mapping[labelx] == 'trained':
                    trained_scores.extend(cv_scores.tolist())                    
                else:
                    untrained_scores.extend(cv_scores.tolist())                   
                    
    accuracy_trained = np.mean(trained_scores)
    accuracy_untrained = np.mean(untrained_scores)
    accuracy_trained_untrained = np.mean(trained_untrained_scores) 

    return(accuracy_trained, accuracy_untrained, accuracy_trained_untrained)   


def get_cov_ledoit(residual_file_list, mymasker, constant_columns):

    d = np.sum(~constant_columns)
    cov_ledoit = np.zeros(d, d)
    
    for residual_file in residual_file_list:
        residuals = mymasker.fit_transform(residual_file)[:, ~constant_columns]
        
        cov_ledoit_partial, _ = ledoit_wolf(residuals)    
        cov_ledoit = cov_ledoit + cov_ledoit_partial
        
    cov_ledoit = cov_ledoit/len(residual_file_list)
    Uc, Dc, Vhc = svd(cov_ledoit, full_matrices = False)
    cov_ledoit_sqrt = np.dot(Uc * np.sqrt(Dc), Vhc)
    
    return (cov_ledoit_sqrt)

def prewhiten(betas, cov_ledoit_sqrt):
#    cov_ledoit, _ = ledoit_wolf(residuals)    
#    Uc, Dc, Vhc = svd(cov_ledoit, full_matrices = False)
#    cov_ledoit_sqrt = np.dot(Uc * np.sqrt(Dc), Vhc)
    prewhiten_ok = True
    try:
        betas_prewhitened = np.dot(betas, np.linalg.inv(cov_ledoit_sqrt))
    except np.linalg.LinAlgError:
        betas_prewhitened = betas #if there are problems use original
        prewhiten_ok = False
    return(betas_prewhitened, prewhiten_ok)

def process(label, subject, session, subjects_dir, analysis_dir, labels_dir, 
            sequences, effects_file, roi_data_file, do_prewhitening = None, 
            overwrite_extract = False,
            permutate = False):
            
    hemi = 'rh' if label.startswith('R_') else 'lh'
    mylabel = project_label(label, hemi, subject, 
                                    subjects_dir, 
                                    analysis_dir, 
                                    labels_dir, 
                                    overwrite = overwrite_extract)
    
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

    sequences.loc[:, 'iscorrect'] = sequences.accuracy >= MINACC
    # to accept the data of a run, there have to be at least 2 correct instances
    # of each sequence type in it
    valid_types = sequences.groupby(['run', 'seq_type']).iscorrect.sum()
    valid_runs = (valid_types >= MINCORRECT).groupby('run').sum() == NTYPES
    valid_runs = valid_runs.index[valid_runs]
    valid = np.logical_and(sequences.run.isin(valid_runs),
                           sequences.iscorrect)
    nvalid = np.sum(valid)
    nn = len(valid)
    nruns = len(valid_runs)

    # remove constant columns
    constant_columns = np.all(effects[1:] == effects[:-1], axis=0)
    
    effects_orig = effects.copy()
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

            
    effects_orig = effects_orig*0
    effects_orig[:, ~constant_columns] = effects
    
    # write everything out
    colnames = ['subject',
                'hemi',
                'label',
                'session',
                'run',
                'true_sequence',
                'seq_train',
                'seq_type',
                'iscorrect', 
                'valid'] + \
        ['feature_%0.4d'%(x) for x in range(effects_orig.shape[1])]
                
    alldata = pd.concat([pd.Series(nn*[subject]), 
                         pd.Series(nn*[hemi]), 
                         pd.Series(nn*[label]), 
                         pd.Series(nn*[session]), 
                         sequences.run, 
                         sequences.true_sequence, 
                         sequences.seq_train,
                         sequences.seq_type, 
                         sequences.iscorrect,
                         valid,
                         pd.DataFrame(effects_orig)], axis = 1)
    
    alldata.columns = colnames
    if roi_data_file != None:
        alldata.to_csv('%s_%s_%s.csv'%(roi_data_file, hemi, label),
                       index = False, 
                       float_format = '%.4f')
        
    sequences = sequences.loc[ valid, : ]

    targets = sequences.seq_type.to_numpy()
    seq_train = sequences.seq_train.to_numpy()
    runs = sequences.run.to_numpy()
    effects = effects[valid, :]
    
    targets_perm = np.array(targets, copy = True)
    
    for run in np.unique(runs):
        targets_perm[runs == run] = \
            np.random.permutation(targets[runs == run])    

    if permutate:

        for run in np.unique(runs):
            myperm = np.random.permutation(range(np.sum(runs == run)))
            targets[runs == run] = targets[runs == run][myperm]
            seq_train[runs == run] = seq_train[runs == run][myperm]


    try:
        # print("Computing SVM classification for label %s"%label)    
        # SVM
        clf_acc = fit_clf(effects, targets, runs)
        df = sequences[['seq_type', 'seq_train']].drop_duplicates()

        clf_acc_trained, clf_acc_untrained,  clf_acc_trained_untrained = \
            fit_clf_separate(effects, targets, runs, df)

        clf_acc_perm = fit_clf(effects, targets_perm, runs)

        mean_signal_trained = \
            np.mean(
                np.mean(effects[sequences.seq_train.to_numpy() == 'trained' , :], 
                        axis = 0))
        mean_signal_untrained = \
            np.mean(
                np.mean(effects[sequences.seq_train.to_numpy() == 'untrained' , :], 
                axis = 0))

    except:  
        clf_acc, clf_cm_acc, clf_acc_trained, clf_acc_untrained, \
            clf_acc_trained_untrained, \
            mean_signal_trained, mean_signal_untrained,\
            clf_acc_perm = [np.nan]*7

    try: 
        
        ##############################
 #       if False:
            # prewhiten the data
#            residual_file_list = [os.path.join(analysis_dir, subject, 
#                                         'ses-%d'%session, 
#                                         'run%d'%run,
#                                         'res4d_LSA_%d.nii.gz'%run) \
#                                  for run in np.unique(sequences.run)]
#            cov_ledoit_sqrt = get_cov_ledoit(residual_file_list, nifti_masker, 
#                                             constant_columns)
 #           for run in np.unique(sequences.run):
 #               effects[sequences.run == run, :], prewhiten_ok = prewhiten(
 #                   effects[sequences.run == run, :], cov_ledoit_sqrt)        
        ##############################

        # compute the different indices
        XG = second_moment(sequences)
        theta, resnorm = fit_second_moment(effects, XG)
        T_ind, theta_PCM, G_hat_trained, G_hat_untrained = \
        second_moment_PCM(sequences, effects)

        alpha_trained = np.log(theta[6])
        alpha_untrained = np.log(theta[7])

        alpha_trained_PCM = theta_PCM[5]
        alpha_untrained_PCM = theta_PCM[6]
    except:
        alpha_trained, alpha_untrained, alpha_trained_PCM,\
            alpha_untrained_PCM, resnorm, \
            G_hat_trained, G_hat_untrained = [np.nan]*7
    
    try:
        xnobis, xnobis_same, xnobis_different, xnobis_trained_same, \
            xnobis_trained_different, xnobis_untrained_same, \
            xnobis_untrained_different, xnobis_trained_untrained, \
            labels = compute_distance(effects, targets, df, splits = runs,
                                      dist_type = 'xnobis')
    
    except:
        xnobis_same, xnobis_different, xnobis_trained_same, \
        xnobis_trained_different, xnobis_untrained_same, \
        xnobis_untrained_different, xnobis_trained_untrained = [np.nan]*7

    try:
        cosine, cosine_same, cosine_different, cosine_trained_same, \
            cosine_trained_different, cosine_untrained_same, \
            cosine_untrained_different, cosine_trained_untrained, \
            labels = compute_distance(effects, targets, df, splits = runs,
                                      dist_type = 'cosine')
    
    except:
        cosine_same, cosine_different, cosine_trained_same, \
        cosine_trained_different, cosine_untrained_same, \
        cosine_untrained_different, cosine_trained_untrained = [np.nan]*7

    return((subject, session, hemi, label, resnorm,
            alpha_trained, 
            alpha_untrained, 
            alpha_trained_PCM, 
            alpha_untrained_PCM, 
            G_hat_trained,
            G_hat_untrained,
            clf_acc,
            clf_acc_perm,
            clf_acc_trained, 
            clf_acc_untrained,
            clf_acc_trained_untrained,
            xnobis_same, 
            xnobis_different, 
            xnobis_trained_same, 
            xnobis_trained_different, 
            xnobis_untrained_same, 
            xnobis_untrained_different, 
            xnobis_trained_untrained,
            cosine_same, 
            cosine_different, 
            cosine_trained_same, 
            cosine_trained_different, 
            cosine_untrained_same, 
            cosine_untrained_different, 
            cosine_trained_untrained,
            mean_signal_trained, 
            mean_signal_untrained,
            nvalid, 
            nruns,
            effects.shape[1], 
            constant_columns.shape[0]))

def gather_results(analysis_dir, suffix):
    files = glob(os.path.join(analysis_dir, 'sub-*', 'ses-*', 
                                   'surf', 'roi_scores_%s.csv'%(suffix))
                      )
    results = []
    for file in files:
        results.append(pd.read_csv(file))
        
    results = pd.concat(results)
    return(results)