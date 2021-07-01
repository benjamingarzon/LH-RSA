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
from sklearn.model_selection import cross_val_score, PredefinedSplit, \
    GridSearchCV, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score, make_scorer
from sklearn.covariance import ledoit_wolf
from scipy.optimize import nnls
from scipy.linalg import svd
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_validate, cross_val_predict

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


param_grid = {
    "estimator__svc__C": [np.exp(i) for i in np.arange(-12, -1, 0.5)]#,
#    "estimator__svc__gamma": [np.exp(i) for i in np.arange(-3, 1, 0.5)]
    }  

#inner_cv = KFold(n_splits=10, shuffle=True)
  
balanced_scorer = make_scorer(balanced_accuracy_score)

svc0 = Pipeline([
        ('variancethreshold', VarianceThreshold()),
        ('standardscaler', StandardScaler()),
        ('anova', SelectKBest(f_classif, k = MINFEATURES)),
        ('svc', SVC(kernel='linear'))
    ])

svc_ova = OneVsRestClassifier(svc0)

svc_ovo = OneVsOneClassifier(svc0)

svc_ovo_nl = OneVsOneClassifier(Pipeline([
        ('standardscaler', StandardScaler()),
        ('anova', SelectKBest(f_classif, k = MINFEATURES)),
        ('svc', SVC(kernel='rbf'))
    ]))

balanced_scorer = make_scorer(balanced_accuracy_score)

def fit_svm(X, y, runs, cl = svc_ovo):
    
    # cl_cv = GridSearchCV(estimator = cl, 
    #                       param_grid=param_grid, 
    #                       cv=PredefinedSplit(runs), verbose = 0, 
    #                       scoring = balanced_scorer)
    
    cv_scores = cross_val_score(cl, X, y, cv = PredefinedSplit(runs), 
                                verbose = 0, 
                                scoring = balanced_scorer)
    #cl_cv.fit(X,y)
    #print(cl)
    #print(cl.cv_results_)
    #print(cl.best_params_)
    #print(cv_scores)
    # plt.plot(param_grid['estimator__svc__C'], 
    #          cl_cv.cv_results_['mean_test_score'])
    # plt.show()
    return(cv_scores.mean())


def fit_svm_confusion(X, y, runs, df, cl = svc_ovo):
    mapping = pd.Series(df.seq_train.values, index=df.seq_type) 
 
    # cl_cv = GridSearchCV(estimator = cl, 
    #                   param_grid=param_grid, 
    #                   cv=PredefinedSplit(runs), verbose = 0, 
    #                   scoring = balanced_scorer)
    
    y_pred = cross_val_predict(cl, X, y, cv = PredefinedSplit(runs),
                                verbose = 0)
    
    conf_mat = confusion_matrix(y, y_pred)
    accuracies = conf_mat.diagonal()/np.sum(conf_mat, axis = 1)
    accuracy_trained = np.mean(accuracies[mapping == 'trained'])
    accuracy_untrained = np.mean(accuracies[mapping == 'untrained'])
    
#    cv_results = cross_validate(cl, X, y, cv = PredefinedSplit(runs),
#                                verbose = 1,
#                                scoring = confusion_matrix_scorer)
#    cv_scores = cross_val_score(cl, X, y, cv = PredefinedSplit(runs), 
#                                verbose = 0, 
#                                scoring = balanced_scorer)

    return(accuracies.mean(), accuracy_trained, accuracy_untrained)   
    
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
            sequences, effects_file, roi_data_file, do_prewhitening = False, 
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
    if roi_data_file != None:
        alldata.to_csv(roi_data_file, index = False, float_format = '%.4f')
    allsequences = sequences
    sequences = sequences.loc[ valid, : ]

    targets = sequences.seq_type.array.copy()
    runs = sequences.run
    effects = effects[valid, :]
    # remove constant columns
    constant_columns = np.all(effects[1:] == effects[:-1], axis=0)
    effects = effects[:, ~constant_columns]
    if permutate:
        for run in np.unique(runs):
            targets[runs == run] = np.random.permutation(targets[runs == run])
    # second moment
    # print("Computing second moment matrix for label %s"%label)

    try:
        # print("Computing SVM classification for label %s"%label)    
        # SVM
        svm_acc = fit_svm(effects, targets, runs)
        df = sequences[['seq_type', 'seq_train']].drop_duplicates()
        svm_cm_acc, svm_acc_trained, svm_acc_untrained = \
            fit_svm_confusion(effects, targets, runs, df)

        mean_signal_trained = \
            np.mean(
                np.mean(effects[sequences.seq_train.to_numpy() == 'trained' , :], 
                        axis = 0))
        mean_signal_untrained = \
            np.mean(
                np.mean(effects[sequences.seq_train.to_numpy() == 'untrained' , :], 
                axis = 0))

    except:  
        svm_acc, svm_cm_acc, svm_acc_trained, svm_acc_untrained, \
            mean_signal_trained, mean_signal_untrained = [np.nan]*6

    # prewhiten the data
    if do_prewhitening:
        # prewhiten the data
        for run in np.unique(sequences.run):
            residuals_file = os.path.join(analysis_dir, subject, 
                                     'ses-%d'%session, 
                                     'run%d'%run,
                                     'res4d_LSA_%d.nii.gz'%run)
            residuals = nifti_masker.fit_transform(residuals_file)    
            effects[sequences.run == run, :], prewhiten_ok = prewhiten(
                effects[sequences.run == run, :], residuals[:, ~constant_columns])        
    try: 
        # compute the different indices
        XG = second_moment(sequences)
        theta, resnorm = fit_second_moment(effects, XG)
        T_ind, theta_PCM, xnobis_trained, xnobis_untrained = \
        second_moment_PCM(sequences, effects)

        alpha_trained = np.log(theta[6])
        alpha_untrained = np.log(theta[7])

        alpha_trained_PCM = theta_PCM[5]
        alpha_untrained_PCM = theta_PCM[6]
    except:
        alpha_trained, alpha_untrained, alpha_trained_PCM,\
            alpha_untrained_PCM, xnobis_trained, xnobis_untrained = [np.nan]*6
    
    return((subject, session, hemi, label, resnorm,
            alpha_trained, 
            alpha_untrained, 
            alpha_trained_PCM, 
            alpha_untrained_PCM, 
            xnobis_trained,
            xnobis_untrained,
            svm_acc, 
            svm_cm_acc, 
            svm_acc_trained, 
            svm_acc_untrained,
            mean_signal_trained, 
            mean_signal_untrained,
            effects.shape[1]))

def gather_results(analysis_dir, suffix):
    files = glob(os.path.join(analysis_dir, 'sub-*', 'ses-*', 
                                   'surf', 'roi_scores_%s.csv'%(suffix))
                      )
    results = []
    for file in files:
        results.append(pd.read_csv(file))
        
    results = pd.concat(results)
    return(results)