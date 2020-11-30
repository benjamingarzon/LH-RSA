# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os, sys
from mvpa2.suite import *
from mvpa2.base.hdf5 import h5save

import pandas as pd

if __debug__:
    from mvpa2.base import debug
    debug.active += ["SVS", "SLC"]

if True:
    # Define surface and volume data paths:
    datapath = sys.argv[1] 
    sequences_fn = sys.argv[2]
    hemi = sys.argv[3] 
    radius = float(sys.argv[4])
    runs = np.array([int(x) for x in sys.argv[5].split(' ') ])    

else:    
    WD = '/data/lv0/MotorSkill/'
    datapath = os.path.join(WD, 'fmriprep/analysis/sub-lue1101/ses-1/')
    sequences_fn = os.path.join(WD, 'responses/sub-lue1101/ses-1/sequences.csv')
    hemi = 'rh'
    radius = 10.0
    runs = np.array([int(x) for x in "1 2 3 4 5".split(' ') ])
    
surfpath = os.path.join(datapath, 'surf')

epi_fn1 = os.path.join(datapath, 'effects.nii.gz')
epi_fn2 = os.path.join(datapath, 'derivatives.nii.gz')

#load engine 
mask_fn = os.path.join(datapath, 'mask.nii.gz')

allsequences = pd.read_csv(sequences_fn, sep = ' ')
sequences = allsequences.loc[ allsequences.run.isin(runs), : ]

fds1 = fmri_dataset(samples = epi_fn1,
               targets = sequences.seq_type, 
               chunks = sequences.run,
               mask = mask_fn)

pial_surf_fn = os.path.join(surfpath, "%s.pial.gii"
                                     % (hemi))
white_surf_fn = os.path.join(surfpath, "%s.white.gii"
                                     % (hemi))
mid_surf_fn = os.path.join(surfpath, "%s.midthickness.ds.gii"
                                     % (hemi))

if True:
    qe = disc_surface_queryengine(radius, 
                                  fds1,
                                  white_surf = white_surf_fn, 
                                  pial_surf = pial_surf_fn,
                                  source_surf = mid_surf_fn, 
                                  nproc = 30)
    suffix = ''

else:
    
    fds1 = fmri_dataset(samples = epi_fn1,
                   targets = sequences.seq_type, 
                   chunks = sequences.run,
                   mask = mask_fn)
    
#    fds2 = fmri_dataset(samples = epi_fn2,
#                   targets = sequences.seq_type, 
#                   chunks = sequences.run,
#                   mask = mask_fn)

    fds = fds1
#    fds = fds1.copy()
#    fds.samples = np.dstack((fds1.samples, fds2.samples))

    qe = disc_surface_queryengine(radius, 
                                  fds, 
                                  white_surf_fn, 
                                  pial_surf_fn)
    suffix = '.multi'                    

mask = qe.voxsel.get_mask()
nifti_mask = qe.voxsel.get_nifti_image_mask()
nifti_mask.to_filename(os.path.join(surfpath, '%s.qe.mask.nii.gz'%(hemi)))    

#save it  
h5save(os.path.join(surfpath, '%s-%.1f-qe%s.gzipped.hdf5'%(hemi, radius, suffix)), qe, compression=9)

