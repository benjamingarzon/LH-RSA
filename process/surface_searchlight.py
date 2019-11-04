
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os, sys
from mvpa2.suite import *
from mvpa2.clfs.svm import LinearCSVMC
from mvpa2.clfs.ridge import RidgeReg
from mvpa2.base.hdf5 import h5save, h5load
import pandas as pd
from mvpa2.measures import rsa
import pylab as pl
from scipy.spatial.distance import squareform
from utils import *
from scipy import stats
import numpy as np
if __debug__:
    from mvpa2.base import debug
    debug.active += ["SVS", "SLC"]

# Define surface and volume data paths:
minacc = 1.0

if True:
    datapath = sys.argv[1] 
    sequences_fn = sys.argv[2] 
    label_fn = sys.argv[3]
    hemi = sys.argv[4] # 'rh'
    radius = float(sys.argv[5])
    NPROC = float(sys.argv[6])
    testdir = sys.argv[7]
    runs = np.array([int(x) for x in sys.argv[8].split(' ') ])
else:    
    datapath = '/home/benjamin.garzon/Data/LeftHand/Lund1/fmriprep/analysis/sub-lue1105/ses-1/'
    sequences_fn = '/home/benjamin.garzon/Data/LeftHand/Lund1/responses/sub-lue1105/ses-1/sequences.csv'
    label_fn = '/home/benjamin.garzon/Data/LeftHand/Lund1/fmriprep/freesurfer/sub-lue1105/label/rh.cortex.label'
    hemi = 'rh'
    radius = 15.0
    NPROC = 10
    testdir = 'metrics'
    runs = np.array([int(x) for x in "1 2 3 4 5".split(' ') ])

labels = pd.read_csv(label_fn, sep = '\s+', skiprows = 2, header = None)

if not os.path.exists(os.path.join(datapath, 'surf', testdir)):
    os.makedirs(os.path.join(datapath, 'surf', testdir))

# add an explanation about what the test does
    
surfpath = os.path.join(datapath, 'surf', testdir)
epi_fn1 = os.path.join(datapath, 'effects.nii.gz')
epi_fn2 = os.path.join(datapath, 'derivatives.nii.gz')

#load engine 
qe = h5load(os.path.join(datapath, 'surf', '%s-%.1f-qe.gzipped.hdf5'%(hemi, radius)))
mask = qe.voxsel.get_mask()
nifti_mask = qe.voxsel.get_nifti_image_mask()
nifti_mask.to_filename(os.path.join(surfpath, 'qe.mask.nii.gz'))

allsequences = pd.read_csv(sequences_fn, sep = ' ')
sequences = allsequences.loc[ allsequences.run.isin(runs), : ]

print("Total correct: %d of %d"%(np.sum(sequences.accuracy >= minacc), len(sequences.accuracy)))


fds1 = fmri_dataset(samples = epi_fn1,
                   targets =  sequences.seq_type, 
                   chunks = sequences.run,
                   mask = mask)

fds2 = fmri_dataset(samples = epi_fn2,
                   targets = sequences.seq_type, 
                   chunks = sequences.run,
                   mask = mask)

fds1.samples = stats.zscore(fds1.samples, axis = None)
fds2.samples = stats.zscore(fds2.samples, axis = None)

#pl.hist(fds1.samples[:, ::10].ravel(), 1000, alpha=0.5, label='effects')
#pl.hist(fds2.samples[:, ::10].ravel(), 1000, alpha=0.5, label='derivatives')
#pl.legend(loc='upper right')

fds = fds1.copy()
fds.samples = np.dstack((fds1.samples, fds2.samples))
fds.sa.accuracy = sequences.accuracy

roi_ids = np.intersect1d(labels[0].ravel(), qe.ids)

seq_train = sequences.loc[:, ["seq_type", "seq_train"]].drop_duplicates()
seq_train = dict(zip(seq_train["seq_type"], seq_train["seq_train"]))

# compute within spread
for metric in ['correlation']:#['correlation', 'euclidean']:
    dsm_ic = Pwithin_spread(
             seq_train = seq_train, 
             mapper = flatten_mapper(), 
             NCOMPS = 10, 
             filter_accuracy = True,
             accuracy = fds.sa.accuracy,
             pairwise_metric = 'correlation',
             square=False)

    sl_rsa_ic = Searchlight(dsm_ic, 
        queryengine = qe,
        nproc = NPROC, 
        roi_ids = roi_ids)
    results_rsa = sl_rsa_ic(fds)                        
    
    spread_path_fn = os.path.join(surfpath, '%s.sl_within_spread_%s_%.1f.func.gii' % (hemi, metric, radius))
    map2gifti2(results_rsa, 
               spread_path_fn, 
               vertices = qe.voxsel.source.nvertices)

if True:
    # compute classification accuracy 
    
    classifiers = {
            'svm': LinearCSVMC()
                }
    mycl = 'svm'
    
    cv = CrossValidation(classifiers[mycl], 
                         NFoldPartitioner(),
                         errorfx=lambda p, 
                         t: np.mean(p == t),
                         enable_ca=['stats'])
    
    sl_cl = Searchlight(cv, 
        queryengine = qe,
        nproc = NPROC, 
        postproc = mean_sample(), 
        roi_ids = roi_ids)
    
    fds_acc = fds1[sequences.accuracy >= minacc, :]
    
    results_cl = sl_cl(fds_acc)                        
     
    acc_path_fn = os.path.join(surfpath, 
                               '%s.sl_acc_%s_%.1f.func.gii' % (hemi, mycl, radius))
    
    map2gifti2(results_cl, 
               acc_path_fn, 
               vertices = qe.voxsel.source.nvertices)
    
    # compute spread
    for metric in ['correlation']:#['correlation', 'euclidean']:
        dsm_ic = Pspread(mapper = flatten_mapper(), 
                 NCOMPS = 10, 
                 filter_accuracy = True,
                 accuracy = fds.sa.accuracy,
                 pairwise_metric = 'correlation',
                 square=False)
    
        sl_rsa_ic = Searchlight(dsm_ic, 
            queryengine = qe,
            nproc = NPROC, 
            roi_ids = roi_ids)
    
        results_rsa = sl_rsa_ic(fds)                        
        
        spread_path_fn = os.path.join(surfpath, 
                                              '%s.sl_spread_%s_%.1f.func.gii' % (hemi, metric, radius))
        map2gifti2(results_rsa, 
                   spread_path_fn, 
                   vertices = qe.voxsel.source.nvertices)

if False:
    # make sure we have at least two elements and we can calculate distances
    correct = sequences.loc[ sequences.accuracy >= minacc, ["seq_type", "seq_train", "true_sequence", "run"] ]
    correct = correct.groupby(["seq_type", "seq_train", "true_sequence", "run"]).size().to_frame(name = 'count').reset_index()
    correct = correct[correct["count"] > 1 ]
    
    correct[["seq_type", "seq_train", "true_sequence"]].drop_duplicates().sort_values("seq_type").to_csv(
        os.path.join(surfpath, 'sequence_types.csv'), 
        header= False, 
        index= False, sep = ' ', 
        float_format='%.2f')
    
    correct.sort_values("seq_type").to_csv(
        os.path.join(surfpath, 'sequence_counts.csv'), 
        header= False, 
        index= False, sep = ' ', 
        float_format='%.2f')
    
    # compute within spread
    for metric in ['correlation']:#['correlation', 'euclidean']:
        dsm_ic = Pwithin_spread(mapper = flatten_mapper(), 
                 NCOMPS = 10, 
                 filter_accuracy = True,
                 accuracy = fds.sa.accuracy,
                 pairwise_metric = 'correlation',
                 square=False)
        
        sl_rsa_ic = Searchlight(dsm_ic, 
            queryengine = qe,
            nproc = NPROC, 
            roi_ids = roi_ids)
        results_rsa = sl_rsa_ic(fds)                        
        
        spread_path_fn = os.path.join(surfpath, '%s.sl_within_spread_%s_%.1f.func.gii' % (hemi, metric, radius))
        map2gifti2(results_rsa, 
                   spread_path_fn, 
                   vertices = qe.voxsel.source.nvertices)
