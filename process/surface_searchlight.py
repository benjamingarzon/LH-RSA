# -*- coding: utf-8 -*-
"""
Surface searchlight of different measures
"""
import os, sys
from mvpa2.suite import *
from mvpa2.clfs.svm import LinearCSVMC
from mvpa2.base.hdf5 import h5save, h5load
import pandas as pd
from utils import *
from scipy import stats
import numpy as np
if __debug__:
  from mvpa2.base import debug
debug.active += ["SVS", "SLC"]

# Definitions
metrics = ['correlation'] #['correlation', 'euclidean']:
minacc = 1.0

if False:
  datapath = sys.argv[1] 
  sequences_fn = sys.argv[2] 
  label_fn = sys.argv[3]
  hemi = sys.argv[4] # 'rh'
  radius = float(sys.argv[5])
  NPROC = float(sys.argv[6])
  testdir = sys.argv[7]
  runs = np.array([int(x) for x in sys.argv[8].split(' ') ])
else:    
  WD = '/home/benjamin.garzon/Data/LeftHand/Lund1/'
  datapath = os.path.join(WD, 'fmriprep/analysis/sub-lue1101/ses-1/')
  sequences_fn = os.path.join(WD, 'responses/sub-lue1101/ses-1/sequences.csv')
  label_fn = os.path.join(WD, 'fmriprep', 
                          'freesurfer/sub-lue1101/label/rh.cortex.label')
  label_fn = os.path.join(WD, 'fmriprep', 
                          'analysis/sub-lue1101/label/rh.somatomotor-mask.ds.label')

  hemi = 'rh'
  radius = 15.0
  NPROC = 10
  testdir = 'metrics'
  runs = np.array([int(x) for x in "1 2 3 4 5".split(' ') ])

labels = pd.read_csv(label_fn, sep = '\s+', skiprows = 2, header = None)

if not os.path.exists(os.path.join(datapath, 'surf', testdir)):
  os.makedirs(os.path.join(datapath, 'surf', testdir))

# load parameters
surfpath = os.path.join(datapath, 'surf', testdir)
epi_fn1 = os.path.join(datapath, 'effects.nii.gz')
epi_fn2 = os.path.join(datapath, 'derivatives.nii.gz')

# load engine 
qe = h5load(os.path.join(datapath, 'surf', 
                         '%s-%.1f-qe.gzipped.hdf5'%(hemi, radius)))
mask = qe.voxsel.get_mask()
nifti_mask = qe.voxsel.get_nifti_image_mask()
nifti_mask.to_filename(os.path.join(surfpath, 'qe.mask.nii.gz'))

# load response data
allsequences = pd.read_csv(sequences_fn, sep = ' ')
sequences = allsequences.loc[ allsequences.run.isin(runs), : ]
ncorrect = np.sum(sequences.accuracy >= minacc)
print("Total correct: %d of %d"%(ncorrect, 
                                 len(sequences.accuracy)))

fds1 = fmri_dataset(samples = epi_fn1,
                    targets =  sequences.seq_type, 
                    chunks = sequences.run,
                    mask = mask)

#fds2 = fmri_dataset(samples = epi_fn2,
#                    targets = sequences.seq_type, 
#                    chunks = sequences.run,
#                    mask = mask)

#fds1.samples = stats.zscore(fds1.samples, axis = None)
#fds2.samples = stats.zscore(fds2.samples, axis = None)

#fds = fds1.copy()
fds = fds1
#fds.samples = np.dstack((fds1.samples, fds2.samples))
fds.sa.accuracy = sequences.accuracy

roi_ids = np.intersect1d(labels[0].ravel(), qe.ids)

seq_train = sequences.loc[:, ["seq_type", "seq_train"]].drop_duplicates()
seq_train = dict(zip(seq_train["seq_type"], seq_train["seq_train"]))

for metric in []: #metrics:
  # compute within-sequence spread
  dsm_ic = Pwithin_spread(
    seq_train = seq_train, 
    mapper = flatten_mapper(), 
    NCOMPS = 10, 
    filter_accuracy = True,
    accuracy = fds.sa.accuracy,
    pairwise_metric = metric,
    square=False)

  sl_rsa_ic = Searchlight(dsm_ic, 
                          queryengine = qe,
                          nproc = NPROC, 
                          roi_ids = roi_ids)
  results_rsa = sl_rsa_ic(fds)                        

  spread_path_fn = os.path.join(
          surfpath, 
          '%s.sl_within_spread_%s_%.1f.func.gii' % (hemi, metric, radius))
  map2gifti2(results_rsa, 
             spread_path_fn, 
             vertices = qe.voxsel.source.nvertices)
  
  # compute overall spread
  dsm_ic = Pspread(mapper = flatten_mapper(), 
                   NCOMPS = 10, 
                   filter_accuracy = True,
                   accuracy = fds.sa.accuracy,
                   pairwise_metric = 'correlation',
                   square = False)

  sl_rsa_ic = Searchlight(dsm_ic, 
                          queryengine = qe,
                          nproc = NPROC, 
                          roi_ids = roi_ids)
  
  results_rsa = sl_rsa_ic(fds)                        
  
  spread_path_fn = os.path.join(
          surfpath,                                 
          '%s.sl_spread_%s_%.1f.func.gii' % (hemi, metric, radius))
  
  map2gifti2(results_rsa, 
             spread_path_fn, 
             vertices = qe.voxsel.source.nvertices)

if False:
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
    
# Component model
if True:    
    sequences_sub = sequences.loc[sequences.accuracy >= minacc]
    np.sum(sequences_sub.accuracy == 1)
    G_list = []
    G_list.append(np.ones((ncorrect, ncorrect)))
    Xrun = indicator(sequences_sub.run)
    G_list.append(np.dot(Xrun, Xrun.T)) # Run effect
    Xtrained = 1* (sequences_sub.seq_train.values == 'trained' ).reshape(-1, 1)
    G_list.append(np.dot(Xtrained, Xtrained.T)) # Trained common
    Xuntrained = 1* (sequences_sub.seq_train.values == 'untrained' ).reshape(-1, 1)
    G_list.append(np.dot(Xuntrained, Xuntrained.T)) # Untrained common
    
    Xseq1 = indicator(sequences_sub.seq_type*Xtrained.ravel(), positive = True) 
    G_list.append(np.dot(Xseq1, Xseq1.T)) # Difference between two trained sequences
    Xseq2 = indicator(sequences_sub.seq_type*Xuntrained.ravel(), positive = True)
    G_list.append(np.dot(Xseq2, Xseq2.T)) # Difference between two untrained sequences
    
    # the last two parameters are the ones of interest
    G_list.append(np.diag(Xtrained.ravel())) # Variability of trained sequences
    G_list.append(np.diag(Xuntrained.ravel())) # Variability of untrained sequences
    
    G = np.dstack(G_list)
    XG = np.hstack([x.reshape(-1, 1) for x in G_list])
    
    dsm_ic = Psecond_moment(
    XG = XG, 
    mapper = flatten_mapper(), 
    filter_accuracy = True,
    accuracy = fds.sa.accuracy,
    return_pars = np.arange(8),
    square=False)

    sl_pcm_ic = Searchlight(dsm_ic, 
                              queryengine = qe,
                              nproc = NPROC, 
                              roi_ids = roi_ids)
    
    results_pcm = sl_pcm_ic(fds)                        
    
    results_pcm_aux = results_pcm.copy()
    
    # compute ratio between untrained and trained in log scale
    results_pcm_aux.samples = 
    np.log(results_pcm.samples[7]/results_pcm.samples[6]).reshape(1, -1)
    pcm_path_fn = os.path.join(
              surfpath, 
              '%s.sl_pcm_%.1f.func.gii' % (hemi, radius))
    map2gifti2(results_pcm_aux, 
                 pcm_path_fn, 
                 vertices = qe.voxsel.source.nvertices)