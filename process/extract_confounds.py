#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 13:37:43 2021

@author: xgarzb@GU.GU.SE
"""

import sys
import pandas as pd
import numpy as np
if True:
    file_in = sys.argv[1]
    file_out = sys.argv[2]
    file_type = sys.argv[3]
else:
    file_in = '/home/xgarzb@GU.GU.SE/Data/LeftHand/Lund1/fmriprep/fmriprep//sub-lue1101/ses-1/func/sub-lue1101_ses-1_task-sequence_run-5_desc-confounds_regressors.tsv'
    file_out = '/home/xgarzb@GU.GU.SE/Data/LeftHand/Lund1/fmriprep/analysis/sub-lue1101/ses-1/run2/confounds.csv'
    file_type = 'extended'

confounds = pd.read_csv(file_in, sep = '\t')

if file_type != 'extended':
    mycols = [x for x in confounds.columns if 'rot_' in x or 'trans_' in x]
#    mycols = ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']
else:
    mycols = [x for x in confounds.columns if 'aroma_motion' in x] # in x or 'rot_' in x or 'trans_' in x]
    #mycols = [x for x in confounds.columns if 'aroma_motion' in x] + \
    #    ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']
        
confounds = confounds.loc[:, mycols]
confounds = confounds.fillna(confounds.mean())
print("Size of {} confound matrix: {}".format(file_type, confounds.shape))
confounds.to_csv(file_out, sep = '\t', header=False, index=False )
