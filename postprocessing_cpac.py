#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

 | Title:
     postprocessing of the output from CPAC
     
 | Date:
     2020-09-20
     
 | Author(s):
     Theodor Rumetshofer

 | Description:    
     This script uses Nipype to extract the define confounds from the CPAC 
     output and standardize the time-series. Additionally high-motion volumes 
     were removed and the output image was smoothed.

 | List of functions:
     No user defined functions are used in the program.

 | List of "non standard" modules:
     None

 | Procedure:
     1) define the confounds and the input data into the pipeline
     2) standardize and detrend the data
     3) resample the rsfMRI and the mask to MNI-space
     4) remove high-motion artefacts by scrubbing
     5) smoothing

 | Usage:
     ./postprocessing_cpac.py

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

#############################################################################


#_______________
# load the data \_____________________________________________________________

# POWEREDGE
PATH = '/mnt/work/u207955/hd/SLE/DATA_7T/rsfmri_200909/'



## OUTLIER: define FD and DVARS 
fd_thresh = 0.4 # only FD is used for the master thesis
#dvars_thresh = 1.5

# number of threads
#num_threads = 6

#_____________________________________
# infosource, input and output stream \_______________________________________

# project folder
project_folder = 'ppCPAC'
# project path
project_path = PATH+project_folder+'/output/output/pipeline_ppCPAC_freq-filter_nuisance/'


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
file_found = project_path+'sub-*/functional_brain_mask_to_standard/_scan_rest/sub-*_brain_mask_maths_antswarp.nii.gz'

# search files
subjects_found = np.sort(gl(file_found))
subject_list = [item.split('/functional_brain_mask_to_standard/')[0].split('/')[-1] for item in subjects_found]

# selector list (NS)
selector_list = ['_selector_WM-2mm-M_CSF-2mm-M_tC-5PCT2-PC5_M-SDB_P-2_BP-B0.01-T0.1',
                 '_selector_WM-2mm-M_CSF-2mm-M_tC-5PCT2-PC5_aC-CSF+WM-2mm-PC5_G-M_M-SDB_P-2_BP-B0.01-T0.1'
                 ]


## infosource - a function free node to iterate over the list of subject names
infosource = Node(IdentityInterface(fields=['subject_id','selector']),name="infosource")

infosource.iterables = [('subject_id', subject_list)
                        , ('selector', selector_list)
                        ]
 

## input stream (SelectFiles)
data_pool = {  
            'bold_mni':'{subject_id}/functional_to_standard/_scan_rest/{selector}/bandpassed_demeaned_filtered_antswarp.nii.gz',
            'bold_mask_mni':'{subject_id}/functional_brain_mask_to_standard/_scan_rest/sub-*_brain_mask_maths_antswarp.nii.gz',


            'fd_power':'{subject_id}/frame_wise_displacement_power/_scan_rest/FD.1D',
            'fd_jenkinson':'{subject_id}/frame_wise_displacement_jenkinson/_scan_rest/FD_J.1D',
            
            'confounds_ts':'{subject_id}/functional_nuisance_regressors/_scan_rest/{selector}/nuisance_regressors.1D',
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
def node_extract_confounds(conf_ts, fd_power, fd_jenkinson):

    import json
    import pandas as pd
    import numpy as np
    import os

    # read all confounds
    confounds_raw = pd.read_csv(conf_ts, delimiter='\t', skiprows=2)
    confounds_raw = confounds_raw.rename(columns={'# RotY':'RotY'})
    
    
    # read FD
    fdp = pd.read_csv(fd_power, header=None)
    fdj = pd.read_csv(fd_jenkinson, header=None)                                              
    fd = pd.concat([fdp, fdj], axis=1)                                              
    fd.columns=['framewise_displacement','rmsd']
    fd = fd.replace(np.nan, 0)
                                          
    # replace nans/infs with 0
    confounds = confounds_raw.replace(np.nan, 0)
        
    return confounds, fd  

# create the node
node_extract_confounds = Node(Function(input_names=['conf_ts','fd_power','fd_jenkinson']
                                    , output_names=['confounds','fd']
                                    , function=node_extract_confounds
                                    ),
                                    name='node_extract_confounds')
    

## save the confound information
def node_save_confounds(confounds, fd):
    
    import os    

    # save confounds for NS
    confound_file = os.path.join(os.getcwd(),'confounds.csv')
    confounds.to_csv(confound_file)

    # save fd
    fd_file = os.path.join(os.getcwd(),'fd.csv')
    fd.to_csv(fd_file)

    
    return confound_file, fd_file
    
    
node_save_confounds = Node(Function(input_names=['confounds','fd']
                                        , output_names=['confound_file','fd_file']
                                        , function=node_save_confounds
                                        ),
                                        name="node_save_confounds")   

      
## standardize data
def node_std(in_file, mask_file):
    
    import json
    import nilearn.image as image
    import nibabel as nib
    import os

    # nuisance regression  
    out_img = image.clean_img(in_file, confounds=None, t_r=None
                              , detrend=True, standardize=True 
#                              , low_pass=0.1, high_pass=0.01
                              , low_pass=None, high_pass=None                              
                              , ensure_finite=True
                              , mask_img=mask_file
                             )              
    out_file = os.path.join(os.getcwd(),'filtered_ns.nii.gz')
    nib.save(out_img, out_file)
    
    
#    out_file = in_file
    
    return out_file

# create the node
node_std = Node(Function(input_names=['in_file','mask_file']
                                    , output_names=['out_file']
                                    , function=node_std
                                    ),
                                    name='node_std')




## RESAMPLE rsfMRI to 1mm isovoxel
node_resample = Node(afni.Resample(
                                  resample_mode = 'Cu'
#                                  , voxel_size= (1.0, 1.0, 1.0)
                                  , outputtype='NIFTI_GZ'
                                  , num_threads = -1
                                  ),
                                  name='node_resample')
#node_resample.iterables = [('resample_mode', ['NN','Li','Cu'])]
    
## RESAMPLE MASK to 1mm isovoxel
node_resample_mask = Node(afni.Resample(
                                        resample_mode = 'NN'
                                        , outputtype='NIFTI_GZ'
                                        , num_threads = -1
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
    

    fd_outliers = (fd.framewise_displacement.values > fd_thresh)*1
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
# defining fd and dvars
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
                        (infosource, selectfiles, [('subject_id', 'subject_id')
                                                    , ('selector', 'selector')
                                                    ])
                        
                        # extract confound values
                        , (selectfiles, node_extract_confounds, [('confounds_ts','conf_ts')
                                                                , ('fd_power','fd_power')
                                                                , ('fd_jenkinson','fd_jenkinson')
                                                                ])
                        
                        # save confound files as csv
                        , (node_extract_confounds, node_save_confounds, [('confounds', 'confounds')
                                                                        , ('fd', 'fd')
                                                                        ])
                        , (node_save_confounds, datasink, [('confound_file', 'confounds_ns')
                                                          , ('fd_file','confounds_fd')
                                                          ])
    
                        # standardize signal
                        , (selectfiles, node_std, [('bold_mni', 'in_file')])
                        , (selectfiles, node_std, [('bold_mask_mni', 'mask_file')
                                                        ])                                          
                        , (node_std, datasink, [('out_file','filtered')])   

                        # resample and change dimension of fMRI
                        , (node_std, node_resample, [('out_file', 'in_file')])
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


analysisflow.run('MultiProc', plugin_args={'n_procs': 6})



print('......THE END')
