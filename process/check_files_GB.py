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

# stages: 

# read scan lists

# freesurfer results

# VBM results

# processed fMRI

# processed MVPA

# processed anatomical

# invalid files 

BIDS_DIR = '~/Data/LeftHand/Lund1/data_BIDS/'
FS_DIR = '~/Data/LeftHand/freesurfer/'
SCAN_LIST_DIR = '/Data/LeftHand/Lund1/QC/scan_lists'
QC_DIR = '~/Data/LeftHand/Lund1/QC/'
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

def check_files(subject, mydir, session, filenames, folder, 
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

## analysis

###############################################################################    
# Read scan_lists and aggregate scan lists
###############################################################################    
scans = []
for f in scan_lists:    
    scans.append(pd.read_csv(f))
    
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
df['complete_all'] = 1*np.logical_and(
        np.logical_and(df.complete_anat == len(anat_files), 
                       df.complete_fmap == len(fmap_files)),
                       df.complete_func == len(func_files))

df.to_csv(os.path.join(QC_DIR, 'BIDS.csv'), index = False)
df[ df.complete_all == False ].to_csv(os.path.join(QC_DIR, 
  'BIDS_incomplete.csv'), index = False)

###############################################################################    
# create BIDS table
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
