#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

 | Title:
     creating figures for the manuscript
     
 | Date:
     2021-02-20
     
 | Author(s):
     Theodor Rumetshofer

 | Description:    
     This script generates the figures which were used in the manuscript

 | List of functions:
     No user defined functions are used in the program.

 | List of "non standard" modules:
     None

 | Procedure:
     The following figures were generated:
         1) head motion from the high-motion subject (not used in the manuscript)
         2) barplot of the scrubbed volumes
         3) confound correlation
         4) correlation matrices of the temporal connectivity

 | Usage:
     ./figures.py

"""

from __future__ import (print_function, division, unicode_literals, absolute_import)     
import numpy as np
import scipy as sp
import pandas as pd
import nibabel as nib
import seaborn as sns
from matplotlib import pyplot as plt
from glob import glob as gl


##############################################################################

path_save = '/Master_report/figures/'


#_______________________________________
# head motion subecjt 132 from fmriprep \______________________________________

hm = pd.read_csv('/Volumes/ID1036/rsfMRI_7T/_ppfMRIprep_filtered/_outputdir/confounds_ns/_ppfMRIprep_filtered/_workingdir/_subject_id_sub-132/_list_confounds_mc24_cos_tcc-5_csfm_wmm/confounds.csv', index_col=0)

trans = hm[['trans_x','trans_y','trans_z']]
rot = hm[['rot_x','rot_y','rot_z']]
rot_degree = rot*180/np.pi

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9,4),  sharex=True)
#fig.suptitle('timated head motion', fontsize=15)
plt.subplots_adjust(wspace=0.3, top=0.85)

axes[0].plot(trans.trans_x, label='x', color='red')
axes[0].plot(trans.trans_y, label='y', color='green')
axes[0].plot(trans.trans_z, label='z', color='blue')
axes[0].set_title('Translation')
axes[0].set_xlabel('volume')
axes[0].set_ylabel('translation (mm)')
axes[0].legend()

axes[1].plot(rot_degree.rot_x, label='x', color='red')
axes[1].plot(rot_degree.rot_y, label='y', color='green')
axes[1].plot(rot_degree.rot_z, label='z', color='blue')
axes[1].set_title('Rotation')
axes[1].set_xlabel('volume')
axes[1].set_ylabel('rotation (°)')
axes[1].legend()

fig.savefig(path_save+'headmotion_sub132_fmriprep.png',bbox_inches='tight', dpi=300)
plt.show()





#_________________________________________________
# store the information from the scrubbed volumes \____________________________

fd_thresh = 0.4

## fmriprep
files = np.sort(gl('/Volumes/ID1036/rsfMRI_7T/_ppfMRIprep_filtered/_outputdir/confounds_fd/_ppfMRIprep_filtered/_workingdir/_subject_id_sub-*/_list_confounds_mc24_cos_tcc-5_csfm_wmm/fd.csv'))
sub_id = [item.split('_subject_id_')[1].split('/')[0] for item in files]

fd_f = pd.DataFrame(columns=sub_id, index=range(0,202))
for file, sub in zip(files, sub_id):
    fd = pd.read_csv(file, index_col=0)
    fd_f[sub] = fd
    
# remove sub 132
fd_f = fd_f.drop('sub-132', axis=1)
# rename
#fd_f.columns = ['fmriprep_'+item for item in fd_f.columns]

## CPAC
files = np.sort(gl('/Volumes/ID1036/rsfMRI_7T/_ppCPAC_filtered/_outputdir/confounds_fd/_ppCPAC_filtered/_workingdir/_selector__selector_WM-2mm-M_CSF-2mm-M_tC-5PCT2-PC5_M-SDB_P-2_BP-B0.01-T0.1_subject_id_sub-*_ses-1/fd.csv'))
sub_id = [item.split('_subject_id_')[1].split('_ses')[0] for item in files]

fd_c = pd.DataFrame(columns=sub_id, index=range(0,202))
for file, sub in zip(files, sub_id):
    fd = pd.read_csv(file, index_col=0)
    fd_c[sub] = fd
    
# remove sub 132
fd_c = fd_c.drop('sub-132', axis=1)


scrub_f_all = []
scrub_c_all = []
for sub1, sub2 in zip(fd_f, fd_c):
    if sub1 == sub2:
        
        scrub_f = np.sum(~(fd_f[sub1].values>fd_thresh))
        scrub_c = np.sum(~(fd_c[sub1].values>fd_thresh))
        scrub_f_all.append(scrub_f)
        scrub_c_all.append(scrub_c)
        
        # get the index of removed volumes
        f = [i for i,x in enumerate(fd_f[sub1].values > fd_thresh) if x]
        c = [i for i,x in enumerate(fd_c[sub2].values > fd_thresh) if x]
        
        # calc also the correlation of fd
        r,p = sp.stats.pearsonr(fd_f[sub1].values, fd_c[sub2].values)
        
        
        print(sub1)
        print('fmriprep:',scrub_f, f)
        print('cpac:',scrub_c, c)
        print(round(r,3), p)
        print()
        



## plot
labels = range(1,10)

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, scrub_f_all, width, label='fMRIprep')
rects2 = ax.bar(x + width/2, scrub_c_all, width, label='CPAC')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('number of volumes after scrubbing')
ax.set_xlabel('subjects')
ax.set_ylim([180, 208])
#ax.set_title('Scores by group and gender')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend(loc='upper left')

#ax.bar_label(rects1, padding=3)
#ax.bar_label(rects2, padding=3)


fig.savefig(path_save+'scrubbing.png',bbox_inches='tight', dpi=300)
plt.show()



#______________________
# confouds correlation \______________________________________________________

## fmriprep
f_files = np.sort(gl('/Volumes/ID1036/rsfMRI_7T/_ppfMRIprep_filtered/_outputdir/confounds_ns/_ppfMRIprep_filtered/_workingdir/_subject_id_sub-*/_list_confounds_mc24_cos_tcc-5_csfm_wmm/confounds.csv'))
subs_f = [item.split('_subject_id_')[1].split('/')[0] for item in f_files]


## CPAC
c_files = np.sort(gl('/Volumes/ID1036/rsfMRI_7T/_ppCPAC_filtered/_outputdir/confounds_ns/_ppCPAC_filtered/_workingdir/_selector__selector_WM-2mm-M_CSF-2mm-M_tC-5PCT2-PC5_M-SDB_P-2_BP-B0.01-T0.1_subject_id_sub-*_ses-1/confounds.csv'))
subs_c = [item.split('_subject_id_')[1].split('_ses')[0] for item in c_files]


## create table
cols = ['trans_x','trans_y','trans_z', 'rot_x', 'rot_y', 'rot_z', 'csf', 'white_matter',
        't_comp_cor_00','t_comp_cor_01','t_comp_cor_02', 't_comp_cor_03', 't_comp_cor_04']

df_table = pd.DataFrame(index=subs_f, columns=cols)


# print results
for sub_f, sub_c, f_file, c_file in zip(subs_f, subs_c, f_files, c_files):
    
    if sub_f == sub_c:
        
        f_conf = pd.read_csv(f_file, index_col=0)
        c_conf = pd.read_csv(c_file, index_col=0)
    
        # take the fmriprep - OBS! search tcomp with regex in case the names of the tcomcpr are different
        f_conf_1 = f_conf[['trans_x','trans_y','trans_z', 'rot_x', 'rot_y', 'rot_z', 'csf', 'white_matter']]
        f_conf_2 = f_conf.filter(regex='t_comp_cor_')
        f_all = pd.concat([f_conf_1, f_conf_2], axis=1)
    
        # take the cpac confounds
        # OBS! Y and Z are vice versa with fmriprep
        c_conf_1 = c_conf[['X','Z','Y', 'RotX', 'RotZ', 'RotY', 'CerebrospinalFluidMean0', 'WhiteMatterMean0']]
        c_conf_2 = c_conf.filter(regex='tCompCorPC')
        c_all = pd.concat([c_conf_1, c_conf_2], axis=1)
    
        print(' |',sub_f)
        for confound_f, confound_c in zip(f_all, c_all):
            r,p = sp.stats.pearsonr(f_all[confound_f], c_all[confound_c])
            print(confound_f, confound_c, np.round(r,3),p, sep='\t')
            
            # save table
            df_table.loc[sub_f][confound_f] = np.round(r,decimals=3)
            
            
        print()
        
## SAVE in TABLE    
# drop subj 132
df_table = df_table.drop('sub-132', axis=0)     

# calc a tCompCor mean + add it to the datafrane
df_table['t_comp_cor_mean'] = df_table.filter(regex='t_comp_cor_').mean(axis=1)

# drop the single t_comp_cor 
df_table = df_table.drop(['t_comp_cor_00','t_comp_cor_01','t_comp_cor_02', 't_comp_cor_03', 't_comp_cor_04'], axis=1)

   
df_table.T.to_excel(path_save+'confounds_corr.xlsx')  
       
        

#__________________
# tFNC correlation \___________________________________________________________

c = '/Volumes/ID1036/rsfMRI_7T/gift210430/output_5/gica_mean_timecourses_ica_s1_.nii'
dc = pd.DataFrame(nib.load(c).get_fdata(), columns=range(1,21))
# tanke only RSN and sort 
dc = dc[[13,2,6,7,20,11,19,3,4,5,12]]
dc = dc.rename(columns={13:'DMN (13)',2:'VIN (2)',6:'VIN (6)', 7:'VIN(7)',
                        20:'lMFN (20)',11:'rMFN (11)',19:'AN (19)',
                        3:'SMN (3)',4:'SMN (4)',5:'SMN (5)',12:'SMN (12)'
                        })
corr_dc = dc.corr()

c_mask = np.tril(np.ones(corr_dc.shape)).astype(np.bool)
hmap = sns.heatmap(corr_dc.round(decimals=2), cmap="Spectral_r", vmax=0.6, vmin=-0.6,
                   annot=True, annot_kws={"size": 8})
hmap.figure.savefig(path_save+'heatmap_cpac.png',bbox_inches='tight', dpi=300)
plt.show()


f = '/Volumes/ID1036/rsfMRI_7T/gift210430/output_5/gica_mean_timecourses_ica_s2_.nii'
df = pd.DataFrame(nib.load(f).get_fdata(), columns=range(1,21))
# tanke only RSN and sort 
df = df[[13,2,6,7,20,11,19,3,4,5,12]]
df = df.rename(columns={13:'DMN (13)',2:'VIN (2)',6:'VIN (6)', 7:'VIN(7)',
                        20:'lMFN (20)',11:'rMFN (11)',19:'AN (19)',
                        3:'SMN (3)',4:'SMN (4)',5:'SMN (5)',12:'SMN (12)'
                        })
corr_df = df.corr()
f_mask = np.tril(np.ones(corr_df.shape)).astype(np.bool)
hmap = sns.heatmap(corr_df.round(decimals=2), cmap="Spectral_r", vmax=0.6, vmin=-0.6,
                   annot=True, annot_kws={"size": 8})
hmap.figure.savefig(path_save+'heatmap_fmriprep.png',bbox_inches='tight', dpi=300)
plt.show()



# Difference

corr_diff = corr_df - corr_dc
mask_diff = np.tril(np.ones(corr_diff.shape)).astype(np.bool)
hmap = sns.heatmap(corr_diff.round(decimals=2), cmap="Spectral_r", vmax=0.6, vmin=-0.6,
                   annot=True, annot_kws={"size": 8})
hmap.figure.savefig(path_save+'heatmap_corr_diff.png',bbox_inches='tight', dpi=300)
plt.show()












print('......THE END')
