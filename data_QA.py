#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

 | Title:
     Quality assessment of the raw data
     
 | Date:
     2020-05-20
     
 | Author(s):
     Theodor Rumetshofer

 | Description:    
     The script generates the follwing images for the QA:
         a) T1
         b) FLAIR
         c) fieldmap
         d) mangitude
         e) rsfMRI (mean image)
      

 | List of functions:
     No user defined functions are used in the program.

 | List of "non standard" modules:
     None

 | Procedure:
     1) Paths of the MRI files in the generated BIDS folder were saved in a dataframe.
     2) for each subjects several QA-images in lightbox format were saved as PNG
     3) all images were linked to a HTML-file for a easier quality assessment

 | Usage:
     ./data_QA.py

"""

from __future__ import (print_function, division, unicode_literals, absolute_import)     
import os
import numpy as np
import scipy as sp
import pandas as pd
import nibabel as nib
from matplotlib import pyplot as plt
from nilearn import image
from glob import glob as gl


###############################################################################


#_______________
# load the data \_____________________________________________________________

# PATH
PATH = '/Users/theo/_analysis/'
img_path = PATH+'SLE/DATA_7T/rsfmri/'

experiment_dir = img_path+'Nifti/'

# create output directory
folder_save = 'QA/0_raw_own/'
path_save = img_path+folder_save 
if not os.path.exists(path_save):
    new_folder = os.makedirs(path_save)   


# T1
df_t1 = experiment_dir+'sub-*/anat/sub-*_T1w.nii.gz'
df_t1 = gl(df_t1)
df_t1 = np.sort(df_t1)
subject_list = [item.split('/Nifti/')[1].split('/')[0] for item in df_t1]
df_t1 = pd.DataFrame(index=subject_list, columns=['t1'], data=df_t1)

# FLAIR
df_flair = experiment_dir+'sub-*/anat/sub-*_FLAIR.nii.gz'
df_flair = gl(df_flair)
df_flair = np.sort(df_flair)
subject_list = [item.split('/Nifti/')[1].split('/')[0] for item in df_flair]
df_flair = pd.DataFrame(index=subject_list, columns=['flair'], data=df_flair)

# rsfMRI
df_rsfmri = experiment_dir+'sub-*/func/sub-*_bold.nii.gz'
df_rsfmri = gl(df_rsfmri)
df_rsfmri = np.sort(df_rsfmri)
subject_list = [item.split('/Nifti/')[1].split('/')[0] for item in df_rsfmri]
df_rsfmri = pd.DataFrame(index=subject_list, columns=['rsfmri'], data=df_rsfmri)

# Fieldmap
df_fieldmap = experiment_dir+'sub-*/fmap/sub-*_phasediff.nii.gz'
df_fieldmap = gl(df_fieldmap)
df_fieldmap = np.sort(df_fieldmap)
subject_list = [item.split('/Nifti/')[1].split('/')[0] for item in df_fieldmap]
df_fieldmap = pd.DataFrame(index=subject_list, columns=['fieldmap'], data=df_fieldmap)

# magnitude
df_magn = experiment_dir+'sub-*/fmap/sub-*_magnitude1.nii.gz'
df_magn = gl(df_magn)
df_magn = np.sort(df_magn)
subject_list = [item.split('/Nifti/')[1].split('/')[0] for item in df_magn]
df_magn = pd.DataFrame(index=subject_list, columns=['magnitude'], data=df_magn)


# put together
df_all = pd.concat([df_t1, df_rsfmri, df_flair, df_fieldmap, df_magn], axis=1, sort=True)

# insert number
df_all.insert(0, value=range(0,df_all.shape[0]), column='num_subjects')

# save
df_all.to_csv(path_save+'list_files.csv')



#________
# run QA \____________________________________________________________________


### run the loop
print('>>> QA is running')
for item in df_all.index:

    df = df_all.loc[item]

    print('  |subject:', df.name)


    #______________________
    # I) T1 - axial slices \__________________________________________________
    
    if isinstance(df['t1'], str): # check if file exists
        title_ = df.name+'_a-t1_axial'
        title_save = path_save+title_
        
        anat = nib.load(df['t1']).get_fdata()
        
        fig, axes = plt.subplots(nrows=4, ncols=6, sharex=True, sharey=True, figsize=(6,6))
        plt.subplots_adjust(hspace=0, wspace=0, top = 0.90, bottom = 0, right = 1, left = 0)
        fig.suptitle(title_, fontsize=15, color='w')
        axes = axes.flatten()
        slices = np.linspace(50, anat.shape[-2]-50, num=len(axes), dtype=int)
        
        for axis, z in zip(axes,slices):
        
            # flip
            axis.imshow(np.flip(sp.ndimage.rotate(anat[:,z,:],180), axis=0), cmap='gray', aspect='auto')    
            
            axis.set_axis_off()
    
        fig.text(-0.02, 0.5, 'LEFT', va='center', rotation='vertical', color='w')
        
        plt.savefig(title_save, dpi=150, facecolor='k', bbox_inches = 'tight', pad_inches = 0)


    #______________________
    # II) FLAIR - axial slices \__________________________________________________
    
    if isinstance(df['flair'], str): # check if file exists
        title_ = df.name+'_b-flair_axial'
        title_save = path_save+title_
        
        flair = nib.load(df['flair']).get_fdata()
        
        fig, axes = plt.subplots(nrows=4, ncols=6, sharex=True, sharey=True, figsize=(6,6))
        plt.subplots_adjust(hspace=0, wspace=0, top = 0.90, bottom = 0, right = 1, left = 0)
        fig.suptitle(title_, fontsize=15, color='w')
        axes = axes.flatten()
        slices = np.linspace(50, flair.shape[-2]-50, num=len(axes), dtype=int)
        
        for axis, z in zip(axes,slices):
        
            # flip
            axis.imshow(np.flip(sp.ndimage.rotate(flair[:,z,:],180), axis=0), cmap='gray', aspect='auto')    
            
            axis.set_axis_off()
    
        fig.text(-0.02, 0.5, 'LEFT', va='center', rotation='vertical', color='w')
        
        plt.savefig(title_save, dpi=150, facecolor='k', bbox_inches = 'tight', pad_inches = 0)


    #_____________________________
    # III) Fieldmap - axial slices \____________________________________________
    
    if isinstance(df['fieldmap'], str): # check if file exists
        
        title_ = df.name+'_c-fieldmap_axial'
        title_save = path_save+title_
        
        fp = nib.load(df['fieldmap']).get_fdata()
        
        fig, axes = plt.subplots(nrows=4, ncols=6, sharex=True, sharey=True, figsize=(6,6))
        plt.subplots_adjust(hspace=0, wspace=0, top = 0.90, bottom = 0, right = 1, left = 0)
        fig.suptitle(title_, fontsize=15, color='w')
        axes = axes.flatten()
        slices = np.linspace(5, fp.shape[-1]-5, num=len(axes), dtype=int)
        
        for axis, z in zip(axes,slices):
        
            # flip
            axis.imshow(np.flip(sp.ndimage.rotate(fp[:,:,z],-90), axis=0), cmap='gray', aspect='auto')    
            
            axis.set_axis_off()
    
        fig.text(-0.02, 0.5, 'LEFT', va='center', rotation='vertical', color='w')
        
        plt.savefig(title_save, dpi=150, facecolor='k', bbox_inches = 'tight', pad_inches = 0)



    #_____________________________
    # IV) Mangitude - axial slices \____________________________________________
    
    if isinstance(df['magnitude'], str): # check if file exists
        
        title_ = df.name+'_d-magnitude_axial'
        title_save = path_save+title_
        
        mag = nib.load(df['magnitude']).get_fdata()
        
        fig, axes = plt.subplots(nrows=4, ncols=6, sharex=True, sharey=True, figsize=(6,6))
        plt.subplots_adjust(hspace=0, wspace=0, top = 0.90, bottom = 0, right = 1, left = 0)
        fig.suptitle(title_, fontsize=15, color='w')
        axes = axes.flatten()
        slices = np.linspace(5, mag.shape[-1]-5, num=len(axes), dtype=int)
        
        for axis, z in zip(axes,slices):
        
            # flip
            axis.imshow(np.flip(sp.ndimage.rotate(mag[:,:,z],-90), axis=0), cmap='gray', aspect='auto')    
            
            axis.set_axis_off()
    
        fig.text(-0.02, 0.5, 'LEFT', va='center', rotation='vertical', color='w')
        
        plt.savefig(title_save, dpi=150, facecolor='k', bbox_inches = 'tight', pad_inches = 0)
    
    #____________
    # V) rsfmri \_____________________________________________________________
    
    if isinstance(df['rsfmri'], str): # check if file exists
        title_ = df.name+'_e-rsfmri_mean_axial'
        title_save = path_save+title_
        
        # load and mean image 
        rs_mean = image.mean_img(df.rsfmri)
        rs_mean = rs_mean.get_fdata()
    #    rs[np.isnan(rs)] = 0
    #    rs_mean = rs.mean(axis=-1)
    
    
        
        fig, axes = plt.subplots(nrows=4, ncols=6, sharex=True, sharey=True, figsize=(6,6))
        plt.subplots_adjust(hspace=0, wspace=0, top = 0.90, bottom = 0, right = 1, left = 0)
        fig.suptitle(title_, fontsize=15, color='w')
        axes = axes.flatten()
        slices = np.linspace(5, rs_mean.shape[-1]-5, num=len(axes), dtype=int)
        
        for axis, z in zip(axes,slices):
        
            # flip
            axis.imshow(np.flip(sp.ndimage.rotate(rs_mean[:,:,z],-90), axis=0), cmap='gray')    
            
            axis.set_axis_off()
    
        fig.text(-0.02, 0.5, 'LEFT', va='center', rotation='vertical', color='w')
        
        plt.savefig(title_save, dpi=150, facecolor='k', bbox_inches = 'tight', pad_inches = 0)
        
    
        
    
    plt.close('all')        
        
#_________________________
# IV) create an html-file \____________________________________________________

 
from dominate import document
from dominate.tags import *

# all pictures
photos = np.sort(gl(path_save+'*.png'))

with document(title='Photos') as doc:
    h1('Photos')
    for path in photos:
        div(img(src=path), _class='photo')

with open(path_save+'index_all.html', 'w') as f:
    f.write(doc.render())    
        

# fieldmap aonly
photos = np.sort(gl(path_save+'*fieldmap_axial.png'))

with document(title='Photos') as doc:
    h1('Photos')
    for path in photos:
        div(img(src=path), _class='photo')

with open(path_save+'index_fieldmaps.html', 'w') as f:
    f.write(doc.render())    


# magntidue aonly
photos = np.sort(gl(path_save+'*magnitude_axial.png'))

with document(title='Photos') as doc:
    h1('Photos')
    for path in photos:
        div(img(src=path), _class='photo')

with open(path_save+'index_magnitudes.html', 'w') as f:
    f.write(doc.render())    


print('......THE END')
