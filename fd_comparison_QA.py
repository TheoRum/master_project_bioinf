#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

 | Title:
     Quality assessment of the framewise displacement (FD)
     
 | Date:
     2020-06-30
     
 | Author(s):
     Theodor Rumetshofer

 | Description:    
     This script extracts and compares the FD extracted from the 10 subjects 
     preprocessed with fMRIprep and CPAC based on criteria by Parkes et al 2018.

 | List of functions:
     No user defined functions are used in the program.

 | List of "non standard" modules:
     None

 | Procedure:
     1) load the confound files for each subject and for each pipeline
     2) extract the FD
     3) calculate a correlation values for each subject between the pipelines
     4) applying the thresholds and rules from Parkes et al 2018
     5) outputs the values and if a subject is fulfilling one or more exclusion
        criteria.

 | Usage:
     ./fd_comparison_QA.py

"""

from __future__ import (print_function, division, unicode_literals, absolute_import)     
             
import numpy as np
import scipy as sp
import pandas as pd
from glob import glob as gl
import json


##############################################################################


#_______________
# load the data \_____________________________________________________________

#PATH = '/Users/theo/_mnt/'
PATH = '/Volumes/ID1036/rsfMRI_7T/' 

# fmriprep
fmriprep_path = PATH+'pp210223/derivatives/fmriprep/'

# CPAC
cpac_path = PATH+'pp210224/output/output/pipeline_pp210224_freq-filter_nuisance/'


#___________________
# define thresholds \_________________________________________________________

# remove complete subjects if: see Parkes et al 2018 (i-iv)

# define FD and DVARS (fmriprep default FD (0.5 mm) and DVARS (1.5))
fd_thresh = 0.4

#mean_fd = 0.5 # Olof
mean_fd = 0.25 # Parkes et al 2018

max_fd = 3 # Olof
#max_fd = 5 # Parkes et al 2018

# max number of subratheshold FDs not more than .. %
#perc_above = 0.2 # Parkes et al 2018

# min length of the data
#min_length = 4 # Parkes et a 2018 
min_length = 5 # Olof

perc_outlier = 0.25 # Parkes et al 2018


print('=> FD for scrubbing:', str(fd_thresh))
print('=> exclusion criteria:')
print(' mean FD:', str(mean_fd))
print(' max FD:', str(max_fd))
print(' length:', str(min_length),'min')
print(' volumes marked as outlier:', str(perc_outlier),'%')



#__________
# fmriprep \___________________________________________________________________

# confounds timesereies
conf_fmriprep = fmriprep_path+'sub-*/func/sub-*_task-rest_desc-confounds_timeseries.tsv'
conf_fmriprep = list(np.sort(gl(conf_fmriprep)))
subjs = [item.split('sub-')[1].split('/')[0] for item in conf_fmriprep] 
conf_fmriprep = pd.DataFrame(conf_fmriprep, columns=['fmriprep'], index=subjs)

# bold TR
bold_tr = fmriprep_path+'sub-*/func/sub-*_task-rest_space-T1w_desc-preproc_bold.json'
bold_tr = list(np.sort(gl(bold_tr)))
subjs = [item.split('sub-')[1].split('/')[0] for item in bold_tr] 
bold_tr = pd.DataFrame(bold_tr, columns=['tr'], index=subjs)


#_______
# c-pac \______________________________________________________________________

# fd
fd_cpac = cpac_path+'sub-*/frame_wise_displacement_power/*/FD.1D'
fd_cpac = list(np.sort(gl(fd_cpac)))
subjs = [item.split('sub-')[1].split('_')[0] for item in fd_cpac] 
fd_cpac = pd.DataFrame(fd_cpac, columns=['fd_cpac'], index=subjs)

# rmsd
rmsd_cpac = cpac_path+'sub-*/frame_wise_displacement_jenkinson/*/FD_J.1D'
rmsd_cpac = list(np.sort(gl(rmsd_cpac)))
subjs = [item.split('sub-')[1].split('_')[0] for item in rmsd_cpac] 
rmsd_cpac = pd.DataFrame(rmsd_cpac, columns=['rmsd_cpac'], index=subjs)



#____
# QA \_________________________________________________________________________

df_all = pd.concat([conf_fmriprep, fd_cpac, rmsd_cpac, bold_tr], axis=1)

corr_r = []
corr_p = []

for item in df_all.index:
        
    df = df_all.loc[item]
    
    # Extract bold TR in seconds
    with open(df.tr, 'r') as fp:
        data = fp.read()
    tr = json.loads(data)['RepetitionTime']    
    
    # read fmirprep  
    conf_fmriprep = df.fmriprep
    conf_fmriprep = pd.read_csv(conf_fmriprep, sep="\t")
    fd_fmirprep = conf_fmriprep.framewise_displacement.values
    dvars_fmriprep = conf_fmriprep.dvars.values
    
    # FD cpac  
    fd_cpac = df.fd_cpac
    fd_cpac = pd.read_csv(fd_cpac, sep="\t").values
    fd_cpac = np.array([item for item in fd_cpac for item in item])
    # rmsd
    rmsd_cpac = df.rmsd_cpac
    rmsd_cpac = pd.read_csv(rmsd_cpac, sep="\t").values
    rmsd_cpac = np.array([item for item in rmsd_cpac for item in item])

    print()
    print('------------------')
    print()
    print(item)
    
    
    ## calc the correlation
    r,p = sp.stats.pearsonr(fd_fmirprep[1:], fd_cpac)
    print('correlation:',r)
    corr_r.append(r)
    corr_p.append(p)
    
    for fd, label in zip([fd_fmirprep,fd_cpac],['FMRIPREP','CPAC']):

        exclude_subj = []
        
        outlier = (fd > fd_thresh)*1
#        dvars_outliers = (conf.std_dvars.values > dvars_thresh)*1
#        outlier = ((dvars_outliers + fd_outliers) > 0)*1
        # default
        #outlier_default = conf_raw.filter(regex='motion_outlier').sum(axis=1).values
        
        
        # i) mean FD > mean_fd
        excl1 = np.nanmean(fd)
        if excl1 > mean_fd: 
            exclude_subj.append(1)
        else: exclude_subj.append(0)
    
        
        # ii) any FD > max_fd
        excl2 = np.sum(fd >= max_fd)
        if excl2 > max_fd:
            exclude_subj.append(1)
        else: exclude_subj.append(0)
        
    
        # iii) > min_length
        excl3 = (np.sum(outlier == 0) * tr)/60
        if excl3 < min_length: 
            exclude_subj.append(1)
        else: exclude_subj.append(0)
        
        # iv) > perc_outlier
        excl4 = 1/len(fd)*np.sum(outlier == 1)
        if excl4 > perc_outlier: 
            exclude_subj.append(1)
        else: exclude_subj.append(0)
        
        # v) total number of outliers from total volume
        excl5 = str(np.sum(outlier == 1))+' / '+str(len(fd))
        
        
        # output
        print()
        print(label)
        print(str(exclude_subj[0]),'| mean FD : ',excl1)
        print(str(exclude_subj[1]),'| volumes with FD >',str(max_fd),'mm :',excl2)    
        print(str(exclude_subj[2]),'| time (min) after scrubbing :' ,excl3)
        print(str(exclude_subj[3]),'| % of scrubbed volumes : ',excl4)
        print('  | total number of scrubbed volume : ',excl5)
        if sum(exclude_subj) > 0:
            print(' => EXCLUDED ')
    






print('......THE END')
