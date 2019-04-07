
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
    datapath = sys.argv[1] #'/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/fmriprep/sub-101/func/'
    sequences_fn = sys.argv[2] #'/home/benjamin.garzon/Data/LeftHand/Lundpilot1/responses/sequences.csv'
    label_fn = sys.argv[3]
    hemi = sys.argv[4] # 'rh'
    radius = float(sys.argv[5]) # 'rh'
    NPROC = float(sys.argv[6]) # 'rh'
    testdir = sys.argv[7]
    runs = np.array([int(x) for x in sys.argv[8].split(' ') ])
else:    
    datapath = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/analysis/sub-105/ses-1/'
    sequences_fn = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/responses/sub-105/ses-1/sequences.csv'
    label_fn = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/freesurfer/sub-105/label/rh.cortex.label'
    hemi = 'rh'
    radius = 15.0
    NPROC = 20
    testdir = 'metrics'
    runs = np.array([int(x) for x in "1 2 3 4".split(' ') ])

labels = pd.read_csv(label_fn, sep = '\s+', skiprows = 2, header = None)

if not os.path.exists(os.path.join(datapath, 'surf', testdir)):
    os.makedirs(os.path.join(datapath, 'surf', testdir))
    # add an explanation to that the test does
    
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

#NPERMS = 100
#permutator = AttributePermutator('targets', count=NPERMS)
#distr_est = MCNullDist(permutator, tail='left', enable_ca=['dist_samples'])

roi_ids = np.intersect1d(labels[0].ravel(), qe.ids)

classifiers = {
#        'elasticnet': 
#             RegressionAsClassifier(clf=SKLLearnerAdapter(sklElasticNet(
#                     alpha=0.01, 
#                     copy_X=True, 
#                     fit_intercept=True, 
#                     l1_ratio=0.3,
#                     max_iter=1000, 
#                     normalize=False, 
#                     positive=False, 
#                     precompute=False,
#                     random_state=None, 
#                     selection='cyclic', 
#                     tol=0.0001, 
#                     warm_start=False), 
#            tags=['enet', 'regression', 'linear', 'does_feature_selection'], 
#            space='targets', 
#            descr='skl.ElasticNet()'), 
#            space='targets', 
#            descr='skl.ElasticNet_C()'),    
        'svm': LinearCSVMC()

#            SVM(svm_impl='NU_SVC', 
#                   kernel=LinearLSKernel(), 
#                   weight=[], probability=1, 
#                   weight_label=[])
            }
mycl = 'svm'

cv = CrossValidation(classifiers[mycl], 
                     NFoldPartitioner(),
                     errorfx=lambda p, 
                     t: np.mean(p == t),
                     enable_ca=['stats'])

if True:
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
    
    np.random.shuffle(fds_acc.targets)
    results_cl = sl_cl(fds_acc)                        
     
    acc_path_fn = os.path.join(surfpath, 
                               '%s.sl_acc_%s_shuffle_%.1f.func.gii' % (hemi, mycl, radius))
    
    map2gifti2(results_cl, 
               acc_path_fn, 
               vertices = qe.voxsel.source.nvertices)

    
#    sl_vol = sphere_searchlight(cv, radius=radius, space='voxel_indices',
#                            nproc = NPROC,
    #                        center_ids = mask.nonzero(),
#                            postproc=mean_sample())
    
#    results_vol = sl_vol(fds1)
    
#    acc_path_fn = os.path.join(surfpath, 
#                               '%s.sl_acc_%s_vol_%.1f.nii.gz' % (hemi, mycl, radius))
    
#    niftiresults = map2nifti(results_vol, imghdr=fds1.a.imghdr)
#    niftiresults.to_filename(acc_path_fn)


#else:
    for mycl in classifiers.keys():
        clf = classifiers[mycl]
        
        cv = CrossValidation(clf, 
                             NFoldPartitioner(),
                             errorfx=lambda p, 
                             t: np.mean(p == t),
                             enable_ca=['stats'])
        
        
        dsm_cl = PClassifier(
                 classifier = cv, 
                 mapper = flatten_mapper(), 
                 NCOMPS = 10, 
                 filter_accuracy = True,
                 accuracy = fds.sa.accuracy)    
        
        sl_cl = Searchlight(dsm_cl, 
            queryengine = qe,
            nproc = NPROC, 
            postproc = mean_sample(), 
            roi_ids = roi_ids)
        
        results_cl = sl_cl(fds)                        
         
        acc_path_fn = os.path.join(surfpath, 
                                   '%s.sl_acc_%s_PCA_%.1f.func.gii' % (hemi, mycl, radius))
        
        map2gifti2(results_cl, 
                   acc_path_fn, 
                   vertices = qe.voxsel.source.nvertices)
    
    
    for metric in ['correlation']:#['correlation', 'euclidean']:
        dsm_ic = Pinvcompactness(mapper = flatten_mapper(), 
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
        
        invcompactness_path_fn = os.path.join(surfpath, 
                                              '%s.sl_invcompactness_%s_%.1f.func.gii' % (hemi, metric, radius))
        map2gifti2(results_rsa, 
                   invcompactness_path_fn, 
                   vertices = qe.voxsel.source.nvertices)
