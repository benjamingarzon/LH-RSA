#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 12:32:11 2020

Script to check which data are complete. 

@author: benjamin.garzon
"""
import glob, os
import numpy as np 
import pandas as pd
import nibabel as nib

# stages: 
# read scan lists
# freesurfer results
# VBM results
# processed fMRI
# processed anatomical
# invalid files 

BIDS_DIR = '/mnt/share/MotorSkill/data_BIDS/'
FS_DIR = '/data/lv0/MotorSkill/freesurfer/'
SCAN_LIST_DIR = '/data/lv0/MotorSkill/QC/scan_lists'
QC_DIR = '/data/lv0/MotorSkill/QC/'
ANALYSIS_DIR = '/data/lv0/MotorSkill/fmriprep/analysis'
FMRIPREP_DIR = '/data/lv0/MotorSkill/fmriprep/fmriprep'

skipped_file = '/data/lv0/MotorSkill/QC/skipped.csv'

# check complete data
scan_lists = glob.glob(os.path.join(SCAN_LIST_DIR, 'Scan_list_wave?.csv'))
BIDS_list = [os.path.basename(x) 
for x in glob.glob(os.path.join(BIDS_DIR, 'sub-*'))]
BIDS_list.sort()
FS_list = [os.path.basename(x) 
for x in glob.glob(os.path.join(FS_DIR, 'sub-*long*base'))]
FS_list.sort()

runs = range(1, 6)

###############################################################################    
# define some functions
###############################################################################    

def check_files_2(subject, mydir, session, filenames, folder, 
                prefix = None, 
                add_ses = True):
    
    if add_ses:
        ses_folder = 'ses-%d'%(session)
    else:
        ses_folder = ''
        
    if prefix:
        filenames = [ prefix +  x for x in filenames ]
            
    available = np.array([
        os.path.exists(
                os.path.join(
                        mydir, 
                        subject, 
                        ses_folder, 
                        folder, 
                        filename.format(subject, session))) 
        for filename in filenames ])
    available = 1*available
    
    return(available.tolist())

def check_files(subject, mydir, session, filenames, folder, 
                prefix = None, 
                add_ses = True,
                count = False):
    
    if add_ses:
        ses_folder = 'ses-%d'%(session)
    else:
        ses_folder = ''
        
    if prefix:
        filenames = [ prefix +  x for x in filenames ]

    def count_vols(filename, count):
        if os.path.exists(filename):
            if count:
                img = nib.load(filename)
                return img.shape[3]
            else:
                return 1
        else:
            return 0
        
        
    available = np.array([
        count_vols(
                os.path.join(
                        mydir, 
                        subject, 
                        ses_folder, 
                        folder, 
                        filename.format(subject, session)),
                count
                ) 
        for filename in filenames ])
    
    return(available.tolist())


def check_fs_files(subject, mydir, session, filenames, folder):
        
    available = np.array([
        os.path.exists(
                os.path.join(
                        mydir, 
                        '{}.{}.long.{}.base'.format(subject, session, subject), 
                        'surf', 
                        filename)) 
        for filename in filenames ])
    #available = 1*np.hstack((available, np.all(available)))
    available = 1*available
    return(available.tolist())

def sumpos(x):
    x.values[x < 0] = 0
    return np.sum(x, axis = 1)

def sumval(x, val):
    return np.sum(x == val, axis = 1)

def itertable(table, runs = None):
    for i, row in table.iterrows():
        if runs:
            yield ('sub-' + row.ID, int(row.SESSION), 
                   np.where([row['FMRI%d'%(run)] != 'nan' for run in runs])[0] + 1)
        else:
            yield 'sub-' + row.ID, int(row.SESSION)
           

###############################################################################    
# define different files
###############################################################################    

## BIDS    
anat_files = [
  'MP2RAGE.nii.gz',
  'MP2RAGE.json',
  'T2w.nii.gz',
  'T2w.json']

fmap_files = [
   'fieldmap.nii.gz', 
   'fieldmap.json',
   'magnitude.nii.gz'
]

fmap_run_files = [ 'run-0%d_fieldmap.nii.gz'%(run) for run in runs ] + \
   [ 'run-0%d_fieldmap.json'%(run) for run in runs ] + \
   [ 'run-0%d_magnitude.nii.gz'%(run) for run in runs ]
# second fmap

BIDS_prefix = '{}_ses-{}_'

func_files = [ 'task-sequence_run-0%d_bold.nii.gz'%(run) for run in runs ] + \
[ 'task-sequence_run-0%d_bold.json'%(run) for run in runs ]

func_invalid_files = [ 'task-sequence_run-0%d_bold_invalid.nii.gz'%(run) for run in runs ] + \
[ 'task-sequence_run-0%d_bold_invalid.json'%(run) for run in runs ]

## freesurfer
fs_files = [
   'lh.pial'
   ]
# analysis
volume_files = [ 'run%d/volume/cope1.nii.gz'%(run) for run in runs ] 

surf_files = [ 'run%d/surfR/cope1.func.gii'%(run) for run in runs ] 

effect_files = [ 'effects.nii.gz' ]

fmriprep_surf_files = [ 'task-sequence_run-%d_space-fsaverage6_hemi-R_bold.func.gii'%(run) for run in runs ] 

fmriprep_vol_files = [ 'task-sequence_run-%d_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'%(run) for run in runs ] 

roi_score_files = [ 'surf/roi_scores.csv' ] 

## analysis

###############################################################################    
# Read scan_lists and aggregate scan lists
###############################################################################    
scans = []
for f in scan_lists:    
    scans.append(pd.read_csv(f, engine = 'python'))
    
mycolumns = ['ID', 'SESSION', 'DIR', 'MP2RAGE', 'MP2RAGE_DelRec', 
             'MP2RAGE_DelRecImag', 'MP2RAGE_DelRecReal', 'B0SHIMMED', 
             'B0SHIMMED2', 'SAG_T2', 'FMRI1', 'FMRI2', 'FMRI3', 'FMRI4', 
             'FMRI5']
scans = pd.concat(scans, sort = False)[mycolumns].astype(str).sort_values(
        ['ID', 'SESSION'])
scans[scans == 'x'] = 'nan'
scans[scans == 'xx'] = 'nan'
scans[scans == 'xxx'] = 'nan'
# clean
scans.to_csv(os.path.join(SCAN_LIST_DIR, 'Scan_list_complete.csv'), 
             index = False)
###############################################################################    
# create BIDS table
###############################################################################    
# code 0: not acquired
#      1: everything ok     
#      -1: acquired but marked as invalid 
#      -2: acquired but missing

results = []
for subject, session, valid_runs in itertable(scans, runs):
    result_anat = check_files(
            subject, BIDS_DIR, session, anat_files, 'anat', BIDS_prefix)
    result_fmap_run = check_files(
            subject, BIDS_DIR, session, fmap_run_files, 'fmap', BIDS_prefix)    
    if np.sum(np.array(result_fmap_run) == 0):
        result_fmap = check_files(subject, BIDS_DIR, session, 
                                  fmap_files, 'fmap', BIDS_prefix)
    else:
        
        result_fmap = [1]*len(fmap_files)        
        print("Several fieldmaps", subject, session)
    result_func = check_files(
            subject, BIDS_DIR, session, func_files, 'func', BIDS_prefix)
    result_invalid_func = check_files(
            subject, BIDS_DIR, session, func_invalid_files, 'func', BIDS_prefix)
 
    f_files = [ 'task-sequence_run-0%d_bold.nii.gz'%(run) for run in valid_runs ] + \
              [ 'task-sequence_run-0%d_bold.json'%(run) for run in valid_runs ]
    
    for f in f_files:
        i = func_files.index(f)
        if result_func[i] == 0: #not found
            result_func[i] = -1 if result_invalid_func[i] == 1 else -2
            
    results.append([subject, 
                    session] + 
                    result_anat + 
                    result_fmap + 
                    result_func)
    
df = pd.DataFrame(results, columns = ['Subject', 'Session'] +  anat_files + 
                  fmap_files + func_files)

# add some summary columns

df['complete_anat'] = sumpos(df[anat_files])
df['complete_fmap'] = sumpos(df[fmap_files])
df['complete_func'] = sumpos(df[func_files])
df['invalid_func'] = sumval(df[func_files], -1)
df['missing_func'] = sumval(df[func_files], -2)
df['complete_BIDS'] = 1*np.logical_and(
        np.logical_and(df.complete_anat == len(anat_files), 
                       df.complete_fmap == len(fmap_files)),
                       df.complete_func == len(func_files))
df_BIDS = df
df_BIDS.to_csv(os.path.join(QC_DIR, 'BIDS.csv'), index = False)
df_BIDSincomplete = df[ df.complete_BIDS == 0 ]
df_BIDSincomplete.to_csv(os.path.join(QC_DIR, 
  'BIDS_incomplete.csv'), index = False)

###############################################################################    
# create FS table
###############################################################################    

results = []
for subject, session in itertable(scans):
    result_fs = check_fs_files(
            subject, FS_DIR, session, fs_files, 'surf')
    results.append([subject, 
                    session] + 
                    result_fs)
        
df = pd.DataFrame(results, columns = ['Subject', 'Session'] +  fs_files)

# add some summary columns
df['complete_fs'] = np.sum(df[fs_files], axis = 1)
df['done'] = 1
df['missing'] = 1 - df['complete_fs']

df.to_csv(os.path.join(QC_DIR, 'FS.csv'), index = False)
df[ df.complete_fs == False ].to_csv(os.path.join(QC_DIR, 
  'FS_incomplete.csv'), index = False)
df_sum = df.groupby(['Subject']).sum()
df_sum.to_csv(os.path.join(QC_DIR, 'FS.sum.csv'), index = True)
df_sum[ df_sum.missing > 0 ].to_csv(os.path.join(QC_DIR, 
      'FS.sum.incomplete.csv'), index = True)
print(df_sum[ df_sum.missing > 0])


###############################################################################    
# create fmriprep table
###############################################################################    
results = []
for subject, session, valid_runs in itertable(scans, runs):
    result_fmriprep_surf = check_files(subject, FMRIPREP_DIR, session, fmriprep_surf_files, 'func', BIDS_prefix)
    result_fmriprep_vol = check_files(subject, FMRIPREP_DIR, session, fmriprep_vol_files, 'func', BIDS_prefix)

    results.append([subject, 
                    session] + 
                    result_fmriprep_surf + 
                    result_fmriprep_vol
                    )
        
df = pd.DataFrame(results, columns = ['Subject', 'Session'] +  fmriprep_surf_files + 
                  fmriprep_vol_files)

df['complete_fmriprep_surf'] = sumpos(df[fmriprep_surf_files])
df['complete_fmriprep_vol'] = sumpos(df[fmriprep_vol_files])
df['complete_FMRIPREP'] = 1*np.logical_and(
    df.complete_fmriprep_surf == len(fmriprep_surf_files),  
    df.complete_fmriprep_vol == len(fmriprep_vol_files)
    )
df = pd.merge(df, df_BIDS[['Subject', 'Session', 'complete_func']], 
                     on = ['Subject', 'Session'])
df['Valid_runs'] = df['complete_func']/2
df['missing_surf'] = - df.complete_fmriprep_surf + df['Valid_runs']
df['missing_vol'] = - df.complete_fmriprep_vol + df['Valid_runs']
df.to_csv(os.path.join(QC_DIR, 'FMRIPREP.csv'), index = False)

df_FMRIPREP = df

df = df[['Subject', 
         'Session', 
         'Valid_runs',
         'complete_fmriprep_surf', 
         'complete_fmriprep_vol', 
         'missing_surf',
         'missing_vol',
         'complete_FMRIPREP']]

df_FMRIPREPincomplete = df[ df.complete_FMRIPREP == False ]
#  np.logical_or(
#    df.missing_surf != 0, 
#    df.missing_vol != 0)
#    ]

df_FMRIPREPincomplete.to_csv(os.path.join(QC_DIR, 
  'FMRIPREP_incomplete.csv'), index = False)

###############################################################################    
# create analysis table
###############################################################################    

results = []
for subject, session, valid_runs in itertable(scans, runs):
    result_scores = check_files(subject, ANALYSIS_DIR, session, roi_score_files, 
                                '')
    result_effects = check_files(subject, ANALYSIS_DIR, session, effect_files, 
                                '', count = True)
    result_volume = check_files(subject, ANALYSIS_DIR, session, volume_files, 
                                '')
    result_surf = check_files(subject, ANALYSIS_DIR, session, surf_files, '')

    results.append([subject, 
                    session] + 
                    result_scores + 
                    result_effects + 
                    result_volume + 
                    result_surf)
  
df = pd.DataFrame(results, columns = ['Subject', 'Session'] + roi_score_files + 
                  effect_files + volume_files + surf_files)

# add some summary columns
df['complete_scores'] = sumpos(df[roi_score_files])
df['complete_effects'] = sumpos(df[effect_files])
df['complete_volume'] = sumpos(df[volume_files])
df['complete_surf'] = sumpos(df[surf_files])

df = pd.merge(df, df_FMRIPREP[['Subject', 'Session', 'complete_func',
                           'complete_fmriprep_surf', 'complete_fmriprep_vol']], 
                     on = ['Subject', 'Session'])

df.to_csv(os.path.join(QC_DIR, 'Analysis.csv'), index = False)

skipped_df = pd.read_csv(skipped_file)

df = pd.merge(df, skipped_df, how = 'left', on = ['Subject', 'Session'])
ntrials = 32
df['enough_trials'] = df.valid_runs.apply(np.isnan)*1
df['complete_analysis'] = 1*np.logical_and.reduce((
                        df.enough_trials == 1,
                        df.complete_scores == 1,
                        df.complete_effects == ntrials*df.complete_fmriprep_vol, 
                        df.complete_volume == df.complete_fmriprep_vol,
                        df.complete_surf == df.complete_fmriprep_surf)
    )

df = df[['Subject', 'Session', 'complete_analysis', 
         'complete_scores',
         'complete_effects', 
         'complete_volume', 'complete_surf',
         'complete_fmriprep_surf', 'complete_fmriprep_vol', 'enough_trials']]


df[ df.complete_analysis == False ].to_csv(os.path.join(QC_DIR, 
  'Analysis_incomplete.csv'), index = False)
# OK : 1201, 3106

