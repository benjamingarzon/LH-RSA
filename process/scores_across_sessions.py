#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 14:56:24 2021

@author: xgarzb@GU.GU.SE
"""
import os
import numpy as np
from roi_analysis_funcs import get_across_session_scores, process_scores
analysis_dir = '/data/lv0/MotorSkill/fmriprep/analysis'
suffix = 'mask-cross'
scores_df = get_across_session_scores(os.path.join(analysis_dir, 
                                'roi_data'),
                                      suffix,
                                      num_cores = 20)

scores_df.to_csv(os.path.join(analysis_dir, 
                                'surf', 
                                'across_session_scores_%s.csv'%(suffix)),
                                index = False,
                                float_format = '%.3f')
