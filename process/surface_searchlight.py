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
  WD = '/data/lv0/MotorSkill/'
  datapath = os.path.join(WD, 'fmriprep/analysis/sub-lue1101/ses-3/')
  sequences_fn = os.path.join(WD, 'responses/sub-lue1101/ses-3/sequences.csv')
  label_fn = os.path.join(WD, 'fmriprep', 
                          'freesurfer/sub-lue1101/label/rh.cortex.label')
#  label_fn = os.path.join(WD, 'fmriprep', 
#                         'analysis/sub-lue1101/label/rh.somatomotor-mask.ds.label')

  hemi = 'rh'
  radius = 10.0
  NPROC = 33
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

NVERTICES = qe.voxsel.source.nvertices
NVERTICES = 40962       

mask = qe.voxsel.get_mask()
#nifti_mask = qe.voxsel.get_nifti_image_mask()
#nifti_mask.to_filename(os.path.join(surfpath, 'qe.mask.nii.gz'))

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

# qe.ids - vertices
roi_ids = np.intersect1d(labels[0].ravel(), qe.ids)

seq_train = sequences.loc[:, ["seq_type", "seq_train"]].drop_duplicates()
seq_train = dict(zip(seq_train["seq_type"], seq_train["seq_train"]))

if True:
    # compute classification accuracy 
    classifiers = {
    'svm': LinearCSVMC()
    }
    mycl = 'svm'
    myresults = []
    
    fds_acc = fds1[sequences.accuracy >= minacc, :].copy()

    for myseq in seq_train.keys():
        print(myseq)
        errorfx = lambda p,t: np.sum(np.logical_and(p == t, t == myseq), dtype = float)/ np.sum(t == myseq, dtype = float)
        cv = CrossValidation(classifiers[mycl], 
        NFoldPartitioner(),
        errorfx = errorfx, 
        enable_ca=['stats'])
        
        sl_cl = Searchlight(cv, 
        queryengine = qe,
        nproc = NPROC,
    #    postproc = mean_sample(), 
        roi_ids = roi_ids)
        
        results_cl = sl_cl(fds_acc)                        
        myresults.append(results_cl)

        acc_path_fn = os.path.join(surfpath, 
        '%s.sl_acc_%s_%.1f_seq-%d.func.gii' % (hemi, mycl, radius, myseq))
        
#        map2gifti2(results_cl, 
#        acc_path_fn, 
#        vertices = NVERTICES)
    # average across sequence types    
    results_cl.samples = np.nanmean(np.concatenate([np.expand_dims(x.samples, 2) for x in myresults], 
                            axis = 2), 
        axis = 2)
    # average across runs
    results_cl = results_cl.get_mapped(mean_sample())
#    results_cl = sl_cl(fds)                        
    
    acc_path_fn = os.path.join(surfpath, 
    '%s.sl_acc_%s_%.1f.func.gii' % (hemi, mycl, radius))
    
    map2gifti2(results_cl, 
    acc_path_fn, 
    vertices = NVERTICES)
#    plot_surf(surf_mesh = os.path.join(surfpath, '../rh.inflated.ds.gii'), 
#                               surf_map = acc_path_fn, hemi = 'right',
#                               view = 'lateral' )
    
#sys.exit()
for metric in metrics: #in metrics:
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
             vertices = NVERTICES)
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
  
  #results_rsa = sl_rsa_ic(fds)                        
  
  spread_path_fn = os.path.join(
          surfpath,                                 
          '%s.sl_spread_%s_%.1f.func.gii' % (hemi, metric, radius))
  
  #map2gifti2(results_rsa, 
  #           spread_path_fn, 
  #           vertices = qe.voxsel.source.nvertices)

#vertex = roi_ids[11000]
#roi_indices = qe.voxsel.volgeom.lin2ijk(qe.voxsel[vertex])
#fd_indices = np.concatenate([ np.where([ np.array_equal(x, y) for x in fds.fa.voxel_indices]) for y in roi_indices ]).ravel()
#cv(fds_acc[:, fd_indices])
#cv(fds1[:, fd_indices])
#nifti_mask = qe.voxsel.get_nifti_image_mask([vertex])
#nifti_mask.to_filename(os.path.join(datapath, 'surf', '%s.%d-roi.nii.gz'%(hemi, vertex)))
                        
  
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

    sl_secmom_ic = Searchlight(dsm_ic, 
                              queryengine = qe,
                              nproc = NPROC, 
                              roi_ids = roi_ids)
    
    results_secmom = sl_secmom_ic(fds)                        
    
    results_secmom_aux = results_secmom.copy()
    
    # compute ratio between untrained and trained in log scale
    results_secmom_aux.samples = \
    np.log(results_secmom.samples[7]/results_secmom.samples[6]).reshape(1, -1)
    secmom_path_fn = os.path.join(
              surfpath, 
              '%s.sl_secmom_%.1f.func.gii' % (hemi, radius))
    map2gifti2(results_secmom_aux, 
                 secmom_path_fn, 
                 vertices = NVERTICES)
