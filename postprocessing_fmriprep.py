#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

 | Title:
     postprocessing of the output from fMRIprep
     
 | Date:
     2020-09-20
     
 | Author(s):
     Theodor Rumetshofer

 | Description:    
     This script uses Nipype to extract the define confounds from the fMRIprep 
     output and apply a nuisance regression. Additionally high-motion volumes 
     were removed and the output image was smoothed.

 | List of functions:
     No user defined functions are used in the program.

 | List of "non standard" modules:
     None

 | Procedure:
     1) define the confounds and the input data into the pipeline
     2) extract the confounds from the files, standardize and detrend the data
     3) resample the rsfMRI and the mask to MNI-space
     4) remove high-motion artefacts by scrubbing
     5) smoothing

 | Usage:
     ./postprocessing_fmriprep.py

"""

from __future__ import (print_function, division, unicode_literals, absolute_import)                  
import time
import os
import numpy as np
import nipype
import nipype.interfaces.afni as afni
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.utility import Function, IdentityInterface
import nipype.interfaces.io as nio 
import nipype as nipy
from glob import glob as gl
from os.path import join as opj


##############################################################################


#_______________
# load the data \_____________________________________________________________

# POWEREDGE
PATH = '/mnt/work/u207955/hd/SLE/DATA_7T/rsfmri_200909/'


## OUTLIER: define FD 
fd_thresh = 0.4 # only FD is used for the master thesis


#___________________________________
# confounds for nuisance regression \_________________________________________

# 24 motion parameters
#mc24 = ['trans_x', 'trans_x_derivative1', 'trans_x_power2',
#'trans_x_derivative1_power2', 'trans_y', 'trans_y_derivative1',
#'trans_y_derivative1_power2', 'trans_y_power2', 'trans_z',
#'trans_z_derivative1', 'trans_z_power2', 'trans_z_derivative1_power2',
#'rot_x', 'rot_x_derivative1', 'rot_x_derivative1_power2',
#'rot_x_power2', 'rot_y', 'rot_y_derivative1',
#'rot_y_derivative1_power2', 'rot_y_power2', 'rot_z',
#'rot_z_derivative1', 'rot_z_power2', 'rot_z_derivative1_power2']

# cosine (detrend and filtering)
#cos = ['cosine00', 'cosine01', 'cosine02', 'cosine03', 'cosine04', 'cosine05']

# WM mask
#wmm = ['white_matter']
#wmmd1 = ['white_matter','white_matter_derivative1']
#wmmd1p = ['white_matter','white_matter_derivative1','white_matter_power2']
#wmmd2p = ['white_matter','white_matter_derivative1','white_matter_power2','white_matter_derivative1_power2']

# CSF mask
#csfm = ['csf']
#csfmd1 = ['csf','csf_derivative1']
#csfmd1p = ['csf','csf_derivative1','csf_power2']
#csfmd2p = ['csf','csf_derivative1','csf_power2','csf_derivative1_power2']

# global signal
#gs = ['global_signal']
#gsd1 = ['global_signal','global_signal_derivative1']
#gsd1p = ['global_signal','global_signal_derivative1','global_signal_power2']
#gsd2p = ['global_signal','global_signal_derivative1','global_signal_power2','global_signal_derivative1_power2']

# COMPCOR - choose x highest components 
# tcompcor 
#tcc-5x
# acompCor CSF
#acccsf-5
# acompCor WM
#accwm-5
# acompCor combined
#acccomb-5

# AROMA - choose all compojnents
#aroma


#list_confounds = [
#        ['mc24','cos','tcc_5','csfm','wmm']
##        , ['mc24','cos','tcc_5','wmm','csfm','gs','acc_comb_5']
##        , ['mc24','cos','csfm','wmm','tcc_5']
##        , ['mc24','cos','csfm_d2p','wmm_d2p','tcc_5'] 
##        , ['mc24','cos','csfm','wmm','gs','tcc_5']
##        , ['mc24','cos','acc_comb_5','csfm','wmm','tcc_5']
##        , ['aroma']
#        ]    



list_confounds = [['mc24_cos_tcc-5_csfm_wmm']]




#_____________________________________
# infosource, input and output stream \_______________________________________

# project folder
project_folder = 'ppfMRIprep'
# project path
project_path = PATH+project_folder+'/derivatives/fmriprep/'


# create output directory
datetimestamp = time.strftime('%Y%m%d%H%M')             
# define the new folder
output_path = PATH+'_'+project_folder+'_filtered'
if not os.path.exists(output_path):
    new_folder = os.makedirs(output_path)   


#  working directory and output folder
output_dir = '_outputdir'
working_dir = '_workingdir'



# location of data folder
file_found = project_path+'sub-*/func/sub-*_task-rest_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz'


# search files
subjects_found = np.sort(gl(file_found))
subject_list = [item.split('fmriprep/')[1].split('/')[0] for item in subjects_found]


## infosource - a function free node to iterate over the list of subject names
infosource = Node(IdentityInterface(fields=['subject_id']),name="infosource")

infosource.iterables = [('subject_id', subject_list)]
 

## input stream (SelectFiles)
data_pool = {
            'bold_mni':'{subject_id}/func/{subject_id}_task-rest_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz',
            'bold_mask_mni':'{subject_id}/func/{subject_id}_task-rest_space-MNI152NLin6Asym_desc-brain_mask.nii.gz',
            'bold_desc':'{subject_id}/func/{subject_id}_task-rest_space-T1w_desc-preproc_bold.json',     
            'confounds_ts':'{subject_id}/func/{subject_id}_task-rest_desc-confounds_timeseries.tsv',
            'confounds_desc':'{subject_id}/func/{subject_id}_task-rest_desc-confounds_timeseries.json',
            }

selectfiles = Node(nipy.SelectFiles(data_pool,
                                    base_directory=project_path),
                        name="selectfiles")

    
# input stream for reference files
data_ref = {'mni_brain': 'MNI_templates/tpl-MNI152NLin2009cAsym_res-01_desc-brain_T1w.nii.gz'
            }
selectrefs = Node(nipy.SelectFiles(data_ref,
                                   base_directory=PATH
                                   ), name='selectrefs')    
    
    
    
    
## output stream (datasink)
datasink = Node(nio.DataSink(base_directory=output_path
                              , container=output_dir),
                    name="datasink")

 
    
#______________
# define Nodes \______________________________________________________________
    
## extract confounds
def node_extract_confounds(conf_desc, conf_ts, list_confounds):

    import json
    import pandas as pd
    import numpy as np
    import os
    from nipype.interfaces.io import ExportFile

    # read all confounds
    confounds_raw = pd.read_csv(conf_ts, sep="\t")
    
    # extract FD
    fd = confounds_raw[['framewise_displacement']] # do it this way otherwise getting an error
    # replace nans/infs with 0
    fd = fd.replace(np.nan, 0)
    
    # extract CompCor+AROMA values
    with open(conf_desc, 'r') as fp:
        data = fp.read()
    obj = json.loads(data)
    
    list_conf = ['name','method','mask','singularvalue','varexpl','retained']
    compcor_aroma = pd.DataFrame(columns=list_conf)
                         
    for o in obj:
        if not 'dropped' in o:
            if 'a_comp_cor_' in o:
                test = {'name':o, 'method':obj[o]['Method'], 'mask':obj[o]['Mask'], 
                        'singularvalue':obj[o]['SingularValue'], 'varexpl':obj[o]['VarianceExplained'], 
                        'retained':obj[o]['Retained']}
                compcor_aroma = compcor_aroma.append(test,ignore_index=True)
            elif 't_comp_cor_' in o:
                test = {'name':o, 'method':obj[o]['Method'], 
                        'singularvalue':obj[o]['SingularValue'], 'varexpl':obj[o]['VarianceExplained'], 
                        'retained':obj[o]['Retained']}
                compcor_aroma = compcor_aroma.append(test,ignore_index=True)
            elif 'aroma_motion_' in o:
                test = {'name':o, 'method':obj[o]['Method']['Name'],
                        'varexpl':obj[o]['TotalVarianceExplained'], 'retained':obj[o]['MotionNoise']}
                compcor_aroma = compcor_aroma.append(test,ignore_index=True)
    # sort df
    compcor_aroma = compcor_aroma.sort_values(by=['method','mask','retained','varexpl'], ascending=False)
    # reindex
    compcor_aroma = compcor_aroma.reset_index(drop=True) 
    
    
    
    new_list_confounds = []

    for item in list_confounds.split('_'):
        
        # Motion
        if 'mc24' in item:
            conf_mc = confounds_raw.filter(regex=r'(trans|rot)').keys().tolist()
            new_list_confounds.append(conf_mc)        
        
        # scanner drift and filtering
        elif 'cos' in item:
            conf_cos = confounds_raw.filter(regex='cosine').keys().tolist()
            new_list_confounds.append(conf_cos)           
        
        # CSF mask
        elif 'csfm' in item:
            if item == 'csfm': conf_csf = ['csf']
            elif item == 'csfmd1': conf_csf = ['csf','csf_derivative1']
            elif item == 'csfmd1p': conf_csf = ['csf','csf_derivative1','csf_power2']
            elif item == 'csfmd2p': conf_csf = ['csf','csf_derivative1','csf_power2','csf_derivative1_power2']
            new_list_confounds.append(conf_csf)              

        # WM mask
        elif 'wmm' in item:
            if item == 'wmm': conf_wm = ['white_matter']
            elif item == 'wmmd1': conf_wm = ['white_matter','white_matter_derivative1']
            elif item == 'wmmd1p': conf_wm = ['white_matter','white_matter_derivative1','white_matter_power2']
            elif item == 'wmmd2p': conf_wm = ['white_matter','white_matter_derivative1','white_matter_power2','white_matter_derivative1_power2']
            new_list_confounds.append(conf_wm) 
 
        # Global signal
        elif 'gs' in item:
            if item == 'gs': conf_gs = ['global_signal']
            elif item == 'gsd1': conf_gs = ['global_signal','global_signal_derivative1']
            elif item == 'gsd1p': conf_gs = ['global_signal','global_signal_derivative1','global_signal_power2']
            elif item == 'gsd2p': conf_gs = ['global_signal','global_signal_derivative1','global_signal_power2','global_signal_derivative1_power2']
            new_list_confounds.append(conf_gs)            
        
        # COMPCOR
        elif 'tcc' in item:
            num_comp = int(item.split('-')[-1])
            conf_tcompcor = compcor_aroma.query('method=="tCompCor"').name[0:num_comp].tolist()
            new_list_confounds.append(conf_tcompcor)
        elif 'acccsf' in item:
            num_comp = int(item.split('-')[-1])
            conf_acompcor_csf = compcor_aroma.query('method=="aCompCor" and mask=="CSF"').name[0:num_comp].tolist()
            new_list_confounds.append(conf_acompcor_csf)
        elif 'accwm' in item:
            num_comp = int(item.split('-')[-1])
            conf_acompcor_wm = compcor_aroma.query('method=="aCompCor" and mask=="WM"').name[0:num_comp].tolist()
            new_list_confounds.append(conf_acompcor_wm)
        elif 'acccomb' in item:
            num_comp = int(item.split('-')[-1])
            conf_acompcor_combined = compcor_aroma.query('method=="aCompCor" and mask=="combined"').name[0:num_comp].tolist()
            new_list_confounds.append(conf_acompcor_combined)
        # AROMA
        elif 'aroma' in item:
            conf_aroma = compcor_aroma.query('method=="ICA-AROMA" and retained==True').name.tolist()
            # rename aroma item because names are not correct in .json and .tsv file
            conf_aroma_rename = []
            for motion in conf_aroma:
                num=motion.split('_')[-1]
                if int(num) < 10:
                    new_num = '_0'+num
                    new_name = '_'.join(motion.split('_')[0:-1]) + new_num
                    conf_aroma_rename.append(new_name)
                else: conf_aroma_rename.append(motion)
            conf_aroma = conf_aroma_rename
            new_list_confounds.append(conf_aroma)  
    
    
    # remodel
    new_list_confounds = [item for item in new_list_confounds for item in item]    # choose confounds    
    confounds = confounds_raw[new_list_confounds]
    # replace nans/infs with 0
    confounds = confounds.replace(np.nan, 0)
        
    return compcor_aroma, confounds, fd 

# create the node
node_extract_confounds = Node(Function(input_names=['conf_desc','conf_ts','list_confounds']
                                    , output_names=['compcor_aroma','confounds','fd']
                                    , function=node_extract_confounds
                                    ),
                                    name='node_extract_confounds')

node_extract_confounds.iterables = ('list_confounds', list_confounds)   
    

## save the confound information
def node_save_confounds(compcor_aroma, confounds, fd):
    
    import os 
    
    # save compCor and aroma
    compcor_aroma_file = os.path.join(os.getcwd(),'compcor_aroma_values.csv')
    compcor_aroma.to_csv(compcor_aroma_file)

    # save confounds for NS
    confound_file = os.path.join(os.getcwd(),'confounds.csv')
    confounds.to_csv(confound_file)

    # fd
    fd_file = os.path.join(os.getcwd(),'fd.csv')
    fd.to_csv(fd_file)
    
    return compcor_aroma_file, confound_file, fd_file
    
    
node_save_confounds = Node(Function(input_names=['compcor_aroma','confounds', 'fd']
                                        , output_names=['compcor_aroma_file','confound_file','fd_file']
                                        , function=node_save_confounds
                                        ),
                                        name="node_save_confounds")   

      
## Nuisance regresseion
def node_nuisance_regression(in_file, mask_file, bold_desc, confounds):
    
    import json
    import nilearn.image as image
    import nibabel as nib
    import os
    
    # Extract bold TR in seconds
    with open(bold_desc, 'r') as fp:
        data = fp.read()
    bold_tr = json.loads(data)['RepetitionTime']

    # nuisance regression  
    out_img = image.clean_img(in_file, confounds=confounds.values, t_r=bold_tr
                              , detrend=True, standardize=True 
                              , low_pass=0.1, high_pass=0.01
#                              , low_pass=None, high_pass=None                              
                              , ensure_finite=True
                              , mask_img=mask_file
                             )              
    out_file = os.path.join(os.getcwd(),'filtered_ns.nii.gz')
    nib.save(out_img, out_file)   
    
    return out_file

# create the node
node_nuisance_regression = Node(Function(input_names=['in_file','mask_file','bold_desc','confounds']
                                    , output_names=['out_file']
                                    , function=node_nuisance_regression
                                    ),
                                    name='node_nuisance_regression')



## RESAMPLE rsfMRI to 1mm isovoxel
node_resample = Node(afni.Resample(
                                  resample_mode = 'Cu'
                                  , outputtype='NIFTI_GZ'
                                  ),
                                  name='node_resample')

    
## RESAMPLE MASK to 1mm isovoxel
node_resample_mask = Node(afni.Resample(
                                        resample_mode = 'NN'
                                        , outputtype='NIFTI_GZ'
                                        ),
                                        name='node_resample_mask')    
    


## SCRUBBING
def node_scrubbing(in_file, fd, fd_thresh):
    
    import json
    import nilearn.image as image
    import nibabel as nib
    import os
    import pandas as pd
    import numpy as np
    

    fd_outliers = (fd.values > fd_thresh)*1
    outlier = (fd_outliers > 0)*1
    
    # INSTEAD: use own code
    outlier_list = [item for item, i in enumerate(outlier) if i==True]
    out_img = nib.load(in_file).get_fdata()
    out_img = np.delete(out_img, outlier_list, axis=3)
    out_img = nib.Nifti1Image(out_img, affine=nib.load(in_file).affine)
    
    
    # name and svae
    out_file = os.path.join(os.getcwd(),'scrubbing.nii.gz')
    nib.save(out_img, out_file)   
    
    return out_file

# create the node
node_scrubbing = Node(Function(input_names=['in_file','fd','fd_thresh']
                                    , output_names=['out_file']
                                    , function=node_scrubbing
                                    ),
                                    name='node_scrubbing')        
# defining fd
node_scrubbing.inputs.fd_thresh = fd_thresh


## Smoothing MERGE 
node_smoothing_afni = Node(afni.Merge(
                        outputtype = 'NIFTI_GZ'
                        , doall = True
                        )
                        , name='node_smoothing_afni')  
node_smoothing_afni.iterables = [('blurfwhm', [2, 3, 4, 5])]



#_________________________________
# define workflow & connect nodes \___________________________________________

#Create a workflow to connect all those nodes
analysisflow = Workflow(name='analysisflow', base_dir=opj(output_path, working_dir))

# connect all nodes to a workflow and with SelectFiles and DataSink to the workflow
analysisflow.connect([
                        # get the source
                        (infosource, selectfiles, [('subject_id', 'subject_id')])
                        
                        # extract confound values
                        , (selectfiles, node_extract_confounds, [('confounds_desc','conf_desc')])
                        , (selectfiles, node_extract_confounds, [('confounds_ts','conf_ts')])
                        
                        # save confound files as csv
                        , (node_extract_confounds, node_save_confounds, [('confounds', 'confounds')
                                                                        , ('compcor_aroma','compcor_aroma')
                                                                        , ('fd', 'fd')
                                                                        ])
                        , (node_save_confounds, datasink, [('compcor_aroma_file', 'confounds_compcor_aroma')
                                                          , ('confound_file', 'confounds_ns')
                                                          , ('fd_file','confounds_fd')
                                                          ])                                          
                        
    
                        # Nuisance regression
                        , (selectfiles, node_nuisance_regression, [('bold_mni', 'in_file')])
                        , (selectfiles, node_nuisance_regression, [('bold_mask_mni', 'mask_file')
                                                                  , ('bold_desc', 'bold_desc')
                                                                  ])                                          
                        , (node_extract_confounds, node_nuisance_regression, [('confounds', 'confounds')])
                        , (node_nuisance_regression, datasink, [('out_file','filtered')])   


                        # resample and change dimension of fMRI
                        , (node_nuisance_regression, node_resample, [('out_file', 'in_file')])
                        , (selectrefs, node_resample, [('mni_brain','master')])
                        , (node_resample, datasink, [('out_file', 'resampled')]) 


                        # resample and change dimension of frmri-mask
                        , (selectfiles, node_resample_mask, [('bold_mask_mni', 'in_file')])
                        , (selectrefs, node_resample_mask, [('mni_brain','master')])
                        , (node_resample_mask, datasink, [('out_file', 'resampled_mask')]) 

    
                        # Scrubbing
                        , (node_resample, node_scrubbing, [('out_file', 'in_file')])
                        , (node_extract_confounds, node_scrubbing, [('fd', 'fd')])
                        , (node_scrubbing, datasink, [('out_file', 'scrubbing')])
                        

                        # Smoothing AFNI
                        , (node_scrubbing, node_smoothing_afni, [('out_file', 'in_files')])
                        , (node_smoothing_afni, datasink, [('out_file', 'smoothing_afni')])    
                    
                        ])


#_________________________
# run workflow / pipeline \___________________________________________________

analysisflow.write_graph(dotfilename='workflow_graph', graph2use='orig', format='png', simple_form=True)

analysisflow.run('MultiProc', plugin_args={'n_procs': 16})



print('......THE END')
