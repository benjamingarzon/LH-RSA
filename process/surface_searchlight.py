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
if __debug__:
    from mvpa2.base import debug
    debug.active += ["SVS", "SLC"]

# Define surface and volume data paths:

if False:
    datapath = sys.argv[1] #'/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/fmriprep/sub-101/func/'
    sequences_fn = sys.argv[2] #'/home/benjamin.garzon/Data/LeftHand/Lundpilot1/responses/sequences.csv'
    hemi = sys.argv[3] # 'rh'
    radius = float(sys.argv[4]) # 'rh'
    NPROC = float(sys.argv[5]) # 'rh'
    
else:    
    datapath = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/analysis/sub-102/'
    sequences_fn = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/responses/sub-102/sequences.csv'
    label_fn = '/home/benjamin.garzon/Data/LeftHand/Lundpilot1/fmriprep/freesurfer/sub-102/label/rh.cortex.label'
    hemi = 'rh'
    radius = 10.0
    labels = pd.read_csv(label_fn, sep = '\s+', skiprows = 2, header = None)
    NPROC = 35
    testdir = 'with_PCA'
    if not os.path.exists(os.path.join(datapath, 'surf', testdir)):
        os.makedirs(os.path.join(datapath, 'surf', testdir))
    # add an explanation
    
surfpath = os.path.join(datapath, 'surf', testdir)
epi_fn1 = os.path.join(datapath, 'effects.nii.gz')
epi_fn2 = os.path.join(datapath, 'derivatives.nii.gz')

#load engine 
qe = h5load(os.path.join(datapath, 'surf', '%s-%.1f-qe.gzipped.hdf5'%(hemi, radius)))
mask = qe.voxsel.get_mask()
nifti_mask = qe.voxsel.get_nifti_image_mask()
nifti_mask.to_filename(os.path.join(surfpath, 'qe.mask.nii.gz'))

sequences = pd.read_csv(sequences_fn, sep = ' ')

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

pl.hist(fds1.samples[:, ::10].ravel(), 1000, alpha=0.5, label='effects')
pl.hist(fds2.samples[:, ::10].ravel(), 1000, alpha=0.5, label='derivatives')
pl.legend(loc='upper right')

fds = fds1.copy()
fds.samples = np.dstack((fds1.samples, fds2.samples))
fds.sa.accuracy = sequences.accuracy

NPERMS = 100
permutator = AttributePermutator('targets', count=NPERMS)
distr_est = MCNullDist(permutator, tail='left', enable_ca=['dist_samples'])


for metric in ['correlation']:#['correlation', 'euclidean']:

#    compact_mapper = FxMapper('samples', 
#                           get_invcompactness_chunk_map, 
#              attrfx='merge', 
#              fxargs = (metric, 
#                       fds.targets, 
#                       fds.chunks))
            
    dsm_multi = PDistMulti(mapper = flatten_mapper(), 
                           pairwise_metric = metric, 
                           square=False)

    dsm_ic = Pinvcompactness(mapper = flatten_mapper(), 
             NCOMPS = 10, 
             filter_accuracy = False,
             pairwise_metric = 'correlation',
             square=False)

    roi_ids = np.intersect1d(labels[0].ravel(), qe.ids)
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
