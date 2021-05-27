#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

 | Title:
     Copying and transformation to BIDS
     
 | Date:
     2020-05-03
     
 | Author(s):
     Theodor Rumetshofer

 | Description:    
     This script copies MRI data (T1, rsfMRI, FLAIR, fieldmap, magnitude) from
     all subjects from a clinical server which belongs to the 017 7T MRI study.
     Additonally, the MRI data were tranformed in the BIDS format by renaming 
     and creating the JSON files. Latter are necessary for the BIDS format and
     the interformation is extracted from the DICOM files.

 | List of functions:
     No user defined functions are used in the program.

 | List of "non standard" modules:
     None

 | Procedure:
     1) search for all the different MRI data on the project folder
     2) save all the paths in a dataframe
     3) save the files in the correct BIDS format
     4) read out the information for the JSON-files from one DICOM-image for 
        each datafile

 | Usage:
     ./data_bids.py

"""
 
import os
import numpy as np
import pandas as pd
from glob import glob as gl
import pydicom
import json


###############################################################################


#______________________________________________________________
# putting all the original and files I analysed in one folder  \_______________
        
## PATHS
PATH = '/Users/data/_analysis/' # mac
img_path = PATH+'SLE/DATA_7T/'


# sink path 
sinkdir = img_path+'rsfmri/Nifti/'
# create folder
if not os.path.exists(sinkdir):
    os.makedirs(sinkdir)  
    
    
# source path
sourcedir = '/Volumes/7T_017_fMRI_SLE/'


# T1
df = np.sort(gl(sourcedir+'7T_017_fMRI_SLE_*/*/NII/*_T13DTFE1mm.nii.gz'))
subject_list = ['7T_017_SLE_'+item.split('_1/NII')[0].split('/')[-2].split('_')[-1] for item in df]
subject_date = [item.split('/NII/')[0].split('/')[-1] for item in df]
serie = [int(item.split('/NII/Serie_')[-1].split('_')[0]) for item in df]
df_t1 = pd.DataFrame(columns=['t1'], data=df)
df_t1.insert(loc=0,value=serie, column='serie_t1')
df_t1.insert(loc=0,value=subject_date, column='date')
df_t1.insert(loc=0,value=subject_list, column='id')

# FLAIR
df = np.sort(gl(sourcedir+'7T_017_fMRI_SLE_*/*/NII/*_3DFLAIR1_ISO_*.nii.gz'))
subject_list = ['7T_017_SLE_'+item.split('_1/NII')[0].split('/')[-2].split('_')[-1] for item in df]
subject_date = [item.split('/NII/')[0].split('/')[-1] for item in df]
serie = [int(item.split('/NII/Serie_')[-1].split('_')[0]) for item in df]
df_flair = pd.DataFrame(columns=['flair'], data=df)
df_flair.insert(loc=0,value=serie, column='serie_flair')
df_flair.insert(loc=0,value=subject_date, column='date')
df_flair.insert(loc=0,value=subject_list, column='id')

# rsfMRI
df = np.sort(gl(sourcedir+'7T_017_fMRI_SLE_*/*/NII/*_WIPFEepi2mm_RS.nii.gz'))
subject_list = ['7T_017_SLE_'+item.split('_1/NII')[0].split('/')[-2].split('_')[-1] for item in df]
subject_date = [item.split('/NII/')[0].split('/')[-1] for item in df]
serie = [int(item.split('/NII/Serie_')[-1].split('_')[0]) for item in df]
df_rs = pd.DataFrame(columns=['rs'], data=df)
df_rs.insert(loc=0,value=serie, column='serie_rs')
df_rs.insert(loc=0,value=subject_date, column='date')
df_rs.insert(loc=0,value=subject_list, column='id')


# Field map
df = np.sort(gl(sourcedir+'7T_017_fMRI_SLE_*/*/NII/*shimmad_32chMTXCLEAR_B0_MAP_B0_UNSPECIFIED.nii.gz'))
subject_list = ['7T_017_SLE_'+item.split('_1/NII')[0].split('/')[-2].split('_')[-1] for item in df]
subject_date = [item.split('/NII/')[0].split('/')[-1] for item in df]
serie = [int(item.split('/NII/')[-1].split('_')[0]) for item in df]
df_fp = pd.DataFrame(columns=['fp'], data=df)
df_fp.insert(loc=0,value=serie, column='serie_fp')
df_fp.insert(loc=0,value=subject_date, column='date')
df_fp.insert(loc=0,value=subject_list, column='id')

# magnitude
df = np.sort(gl(sourcedir+'7T_017_fMRI_SLE_*/*/NII/*shimmad_32chMTXCLEAR_M_FFE_M_FFE.nii.gz'))
subject_list = ['7T_017_SLE_'+item.split('_1/NII')[0].split('/')[-2].split('_')[-1] for item in df]
subject_date = [item.split('/NII/')[0].split('/')[-1] for item in df]
serie = [int(item.split('/NII/')[-1].split('_')[0]) for item in df]
df_mag = pd.DataFrame(columns=['mag'], data=df)
df_mag.insert(loc=0,value=serie, column='serie_mag')
df_mag.insert(loc=0,value=subject_date, column='date')
df_mag.insert(loc=0,value=subject_list, column='id')


### MERGE THE DATA TO THE rsfMRI

##  keep those rs with the newest date
df_rs = df_rs.sort_values(['id','date'])
df_rs = df_rs[~df_rs.duplicated(subset=['id'], keep='last')]

## ADD FIELDMAP
# lets find those fieldmaps that have the same date and the series number is close (max +100) before the rsfMRI series (otherwise it was measured after the rsfMRI)
df_all = df_rs.merge(df_fp, on=['id','date'], how='left')
# get a subset of subjects that are double
subset = df_all[df_all.duplicated('id', keep=False)]
# comare the series within 100
comp_series = subset.serie_rs != subset.serie_fp + 100
remove_subj = subset[comp_series]
# remove the subjects
df_all = df_all.drop(index=remove_subj.index)


## ADD MAGNITUDE
df_all = df_all.merge(df_mag, on=['id','date'], how='left')
# keep those with the same fieldmap serie and magnitude serie
subset = df_all[df_all.duplicated('id', keep=False)]
# comare the series of fieldmap and magnitude
comp_series = subset.serie_fp != subset.serie_mag
remove_subj = subset[comp_series]
# remove the subjects
df_all = df_all.drop(index=remove_subj.index)


## ADD T1w
# keep those the newest T1s with the same date
df_t1 = df_t1[~df_t1.duplicated(subset=['id','date'], keep='last')]
# get all T1s and fieldmaps and magn with the same date as rsfmri image
df_all = df_all.merge(df_t1, on=['id','date'])


# keep the newest FLAIR
df_flair = df_flair.sort_values(['id','date'])
df_flair = df_flair[~df_flair.duplicated(subset=['id','date'], keep='last')]

df_all = df_all.merge(df_flair, on=['id','date'], how='left')

# check for not the same dates and missing values
date_equal = list((df_all.isna().sum(axis=1) == 0)*1)
df_all.insert(loc=2,value=date_equal, column='same_date')


# find the missing flairs
idxs = df_all.query('flair!=flair').index.values
subs = df_all.query('flair!=flair').id.values

for sub, idx in zip(subs, idxs):
    
    # check if FLAIR exists
    if df_flair.query('id==@sub').flair.values.size > 0:
        df_all.set_value(idx, 'flair', df_flair.query('id==@sub').flair.values[0])
        df_all.set_value(idx, 'serie_flair', df_flair.query('id==@sub').serie_flair.values[0])
        
df_all.index = df_all.id.values
df_all = df_all.drop(columns='id')
# insert number
df_all.insert(0, value=range(0,df_all.shape[0]), column='num_subjects')

# save
df_all.to_csv(img_path+'7T_file_list.csv')


### run the loop
print('>>> Copying process is running ...')
# create a dataframe to save the dicom information
df_dcm_rsfmri_all = pd.DataFrame()
df_dcm_t1_all = pd.DataFrame()
df_dcm_flair_all = pd.DataFrame()
df_dcm_fp_all = pd.DataFrame()
df_dcm_mag_all = pd.DataFrame()

for item in df_all.index:

    df = df_all.loc[item]
    print('  |subject:', df.name)
    
    sub_folder = 'sub-'+df.name.split('_')[-1]
    
    ## save 
    dst_path_anat = sinkdir+sub_folder+'/anat/'
    if not os.path.exists(dst_path_anat):
        os.makedirs(dst_path_anat)      
   
    dst_path_func = sinkdir+sub_folder+'/func/'
    if not os.path.exists(dst_path_func):
        os.makedirs(dst_path_func)       
    
    dst_path_fmap = sinkdir+sub_folder+'/fmap/'
    if not os.path.exists(dst_path_fmap):
        os.makedirs(dst_path_fmap)            

   
    # COPY THE FILES + CREATE JSON-FILE
    # delete the dates for the loop
    df = df.drop(['num_subjects','date','same_date','serie_rs','serie_fp','serie_mag','serie_t1','serie_flair'])
    
    ### FLAIR - save dicom header information 
    if isinstance(df.flair, str): # check if file exists
        dcm = gl(df.flair.replace('/NII/','/DCM/').replace('.nii.gz','/*.dcm'))[0]
        hd = pydicom.dcmread(dcm)    
    
        d_flair = {
            'Manufacturer' : hd[0x0008,0x0070].value,
            'ManufacturersModelName' : hd[0x0008,0x1090].value,
            'MagneticFieldStregth' : hd[0x0018,0x0087].value.real,
            'InstitutionName' : hd[0x0008,0x0080].value,
            'PatientsName' : str(hd[0x0010,0x0010].value).replace('_','-'),
            'AcquisitionDate' : hd[0x0008,0x0022].value,
            'MRAcquisitionType' : hd[0x0018, 0x0023].value,
            'SpacingBetweenSlices' : hd[0x0018, 0x0088].value.real,
            'AcquisitionDuration' : hd[0x0018, 0x9073].value.real,
            'RepetitionTime' : hd[0x0018,0x0080].value.real/1000,
            'EchoTime' : hd[0x0018,0x0081].value.real/1000,
            'FlipAngle' : hd[0x0018, 0x1314].value.real,
            'pixdim_1' : hd[0x0028, 0x0030].value[0].real,
            'pixdim_2' : hd[0x0028, 0x0030].value[1].real,
            'pixdim_3' : hd[0x0018, 0x0050].value.real,
            'dim_1' : hd[0x0028, 0x0010].value.real,
            'dim_2' : hd[0x0028, 0x0011].value.real,
            'dim_3' : hd[0x2001, 0x1018].value.real,
            'dim_4' : hd[0x0020, 0x0105].value.real,
            'ParallelReductionFactorInPlane' : hd[0x2005, 0x140f][0][0x0018, 0x9069].value.real,
            'ParallelAcquisitionTechnique' : hd[0x2005, 0x140f][0][0x0018,0x9078].value
            }   
        
        # save in dataframe
        df_dcm_flair = pd.DataFrame.from_dict(data=d_flair, orient='index').T
        df_dcm_flair.index = [df_dcm_flair.PatientsName.values[0].replace('-','_')]
        df_dcm_flair_all = pd.concat([df_dcm_flair_all, df_dcm_flair], sort=False)   
           
   
    ### T1 (for rsfmri) - save dicom header information 
    if isinstance(df.t1, str): # check if file exists
        dcm = gl(df.t1.replace('/NII/','/DCM/').replace('.nii.gz','/*.dcm'))[0]
        hd = pydicom.dcmread(dcm)    
    
        d_t1 = {
            'Manufacturer' : hd[0x0008,0x0070].value,
            'ManufacturersModelName' : hd[0x0008,0x1090].value,
            'MagneticFieldStregth' : hd[0x0018,0x0087].value.real,
            'InstitutionName' : hd[0x0008,0x0080].value,
            'PatientsName' : str(hd[0x0010,0x0010].value).replace('_','-'),
            'AcquisitionDate' : hd[0x0008,0x0022].value,
            'MRAcquisitionType' : hd[0x0018, 0x0023].value,
            'SpacingBetweenSlices' : hd[0x0018, 0x0088].value.real,
            'AcquisitionDuration' : hd[0x0018, 0x9073].value.real,
            'RepetitionTime' : hd[0x0018,0x0080].value.real/1000,
            'EchoTime' : hd[0x0018,0x0081].value.real/1000,
            'FlipAngle' : hd[0x0018, 0x1314].value.real,
            'pixdim_1' : hd[0x0028, 0x0030].value[0].real,
            'pixdim_2' : hd[0x0028, 0x0030].value[1].real,
            'pixdim_3' : hd[0x0018, 0x0050].value.real,
            'dim_1' : hd[0x0028, 0x0010].value.real,
            'dim_2' : hd[0x0028, 0x0011].value.real,
            'dim_3' : hd[0x2001, 0x1018].value.real,
            'dim_4' : hd[0x0020, 0x0105].value.real,
            'ParallelReductionFactorInPlane' : hd[0x2005, 0x140f][0][0x0018, 0x9069].value.real,
            'ParallelAcquisitionTechnique' : hd[0x2005, 0x140f][0][0x0018,0x9078].value
            }   
        
        # save in dataframe
        df_dcm_t1 = pd.DataFrame.from_dict(data=d_t1, orient='index').T
        df_dcm_t1.index = [df_dcm_t1.PatientsName.values[0].replace('-','_')]
        df_dcm_t1_all = pd.concat([df_dcm_t1_all, df_dcm_t1], sort=False)   
                  
    ### rsfmri - save dicom header information 
    if isinstance(df.rs, str): # check if file exists
        dcm = gl(df.rs.replace('/NII/','/DCM/').replace('.nii.gz','/*.dcm'))[0]
        hd = pydicom.dcmread(dcm)    
    
        d_rsfmri = {
            'Manufacturer' : hd[0x0008,0x0070].value,
            'ManufacturersModelName' : hd[0x0008,0x1090].value,
            'MagneticFieldStregth' : hd[0x0018,0x0087].value.real,
            'InstitutionName' : hd[0x0008,0x0080].value,
            'PatientsName' : str(hd[0x0010,0x0010].value).replace('_','-'),
            'AcquisitionDate' : hd[0x0008,0x0022].value,
            'MRAcquisitionType' : hd[0x0018, 0x0023].value,
            'SpacingBetweenSlices' : hd[0x0018, 0x0088].value.real,
            'AcquisitionDuration' : hd[0x0018, 0x9073].value.real,
            'RepetitionTime' : hd[0x0018,0x0080].value.real/1000,
            'EchoTime' : hd[0x0018,0x0081].value.real/1000,
            'FlipAngle' : hd[0x0018, 0x1314].value.real,
            'pixdim_1' : hd[0x0028, 0x0030].value[0].real,
            'pixdim_2' : hd[0x0028, 0x0030].value[1].real,
            'pixdim_3' : hd[0x0018, 0x0050].value.real,
            'dim_1' : hd[0x0028, 0x0010].value.real,
            'dim_2' : hd[0x0028, 0x0011].value.real,
            'dim_3' : hd[0x2001, 0x1018].value.real,
            'dim_4' : hd[0x0020, 0x0105].value.real,
            'ParallelReductionFactorInPlane' : hd[0x2005, 0x140f][0][0x0018, 0x9069].value.real,
            'ParallelAcquisitionTechnique' : hd[0x2005, 0x140f][0][0x0018,0x9078].value,
            'PartialFourier': hd[0x2005, 0x140f][0][0x0018, 0x9081].value,
            'PartialFourierDirection': hd[0x2005, 0x140f][0][0x0018, 0x9036].value,
            # manual etries
            'TaskName' : 'rest',
            # this is "anterior-posterior" but not saved in dicom (see: https://mrtrix.readthedocs.io/en/latest/concepts/pe_scheme.html#fixed-phase-encoding)
            'PhaseEncodingDirection': 'j',     
            # no MB (according to Peter)
            'MultibandAccelerationFactor': 1  
            # remove volumes before pp
#            NumberOfVolumesDiscardedByUser:4
            
            }
        
        ## calc SliceTiming based on the assuption that it is "ascending"
        delta_t = d_rsfmri['RepetitionTime']*1000/(d_rsfmri['dim_3']/1)
        slice_order_odd = list(np.arange(1,d_rsfmri['dim_3']+1,2))
        slice_order_even = list(np.arange(2,d_rsfmri['dim_3']+1,2))
        slice_order = slice_order_odd + slice_order_even
        slice_timing = list(np.multiply(np.arange(0, d_rsfmri['dim_3']),delta_t)/1000)
        df_timing = pd.DataFrame([slice_order, slice_timing]).T
        df_timing = df_timing.sort_values(0)
        
        ## calc the Effective Echo Time and Total ReadoutTime (for the SDC)
        WaterFatShift_ppm = 3.35
        WaterFatShift_ppx = hd[0x2001, 0x1022].value.real
        ImagingFrequency_MHz = hd[0x0018, 0x0084].value.real
        EPIFactor = hd[0x2001, 0x1013].value.real
        EchoTrainLength = EPIFactor + 1
        WaterFatShift_Hz = WaterFatShift_ppm * ImagingFrequency_MHz
        BandWidth_Hz_ppx = WaterFatShift_Hz / WaterFatShift_ppx
        BandWith_tot = BandWidth_Hz_ppx * EchoTrainLength
        EffectiveEchoTime = (1000/(BandWith_tot * d_rsfmri['ParallelReductionFactorInPlane']))/1000 #divide by 1000 to be in seconds [normally 0.3 and 1 ms/acc]
        TotalReadoutTime = EffectiveEchoTime * EchoTrainLength
        
        
        # add to df
        d_rsfmri.update({'SliceTiming': df_timing[1].tolist(),
                         'EffectiveEchoSpacing': EffectiveEchoTime,
                         'TotalReadoutTime': TotalReadoutTime
                         })        

        ## create the JSON file and save it in the func folder
        with open(dst_path_func+sub_folder+'_task-rest_bold.json', 'w') as fp:
            json.dump(d_rsfmri, fp, sort_keys = True, indent = 4,
                   ensure_ascii = False)   
    
        ## save in dataframe
        df_dcm_rsfmri = pd.DataFrame.from_dict(data=d_rsfmri, orient='index').T
        df_dcm_rsfmri.index = [df_dcm_rsfmri.PatientsName.values[0].replace('-','_')]
        df_dcm_rsfmri_all = pd.concat([df_dcm_rsfmri_all, df_dcm_rsfmri], sort=False)   
 
    
           
    ### Magntidue - save dicom header information 
    if isinstance(df.mag, str): # check if file exists
        dcm = gl(df.mag.replace('/NII/','/DCM/').replace('.nii.gz','/*.dcm'))[0]
        hd = pydicom.dcmread(dcm)    
    
        d_mag = {
            'Manufacturer' : hd[0x0008,0x0070].value,
            'ManufacturersModelName' : hd[0x0008,0x1090].value,
            'MagneticFieldStregth' : hd[0x0018,0x0087].value.real,
            'InstitutionName' : hd[0x0008,0x0080].value,
            'PatientsName' : str(hd[0x0010,0x0010].value).replace('_','-'),
            'AcquisitionDate' : hd[0x0008,0x0022].value,
            'MRAcquisitionType' : hd[0x0018, 0x0023].value,
            'SpacingBetweenSlices' : hd[0x0018, 0x0088].value.real,
            'AcquisitionDuration' : hd[0x0018, 0x9073].value.real,
            'RepetitionTime' : hd[0x0018,0x0080].value.real/1000,
            'EchoTime' : hd[0x0018,0x0081].value.real/1000,
            'FlipAngle' : hd[0x0018, 0x1314].value.real,
            'pixdim_1' : hd[0x0028, 0x0030].value[0].real,
            'pixdim_2' : hd[0x0028, 0x0030].value[1].real,
            'pixdim_3' : hd[0x0018, 0x0050].value.real,
            'dim_1' : hd[0x0028, 0x0010].value.real,
            'dim_2' : hd[0x0028, 0x0011].value.real,
            'dim_3' : hd[0x2001, 0x1018].value.real,
            'dim_4' : hd[0x0020, 0x0105].value.real,
            'ParallelReductionFactorInPlane' : hd[0x2005, 0x140f][0][0x0018, 0x9069].value.real,
            'ParallelAcquisitionTechnique' : hd[0x2005, 0x140f][0][0x0018,0x9078].value
            }   
        
        # save in dataframe
        df_dcm_mag = pd.DataFrame.from_dict(data=d_mag, orient='index').T
        df_dcm_mag.index = [df_dcm_mag.PatientsName.values[0].replace('-','_')]
        df_dcm_mag_all = pd.concat([df_dcm_mag_all, df_dcm_mag], sort=False)   
           
            
    
    ### Fieldmap - save dicom header information 
    if isinstance(df.fp, str): # check if file exists
        dcm = gl(df.fp.replace('/NII/','/DCM/').replace('.nii.gz','/*.dcm'))[0]
        hd = pydicom.dcmread(dcm)    
    
        d_fp = {
            'Manufacturer' : hd[0x0008,0x0070].value,
            'ManufacturersModelName' : hd[0x0008,0x1090].value,
            'MagneticFieldStregth' : hd[0x0018,0x0087].value.real,
            'InstitutionName' : hd[0x0008,0x0080].value,
            'PatientsName' : str(hd[0x0010,0x0010].value).replace('_','-'),
            'AcquisitionDate' : hd[0x0008,0x0022].value,
            'MRAcquisitionType' : hd[0x0018, 0x0023].value,
            'SpacingBetweenSlices' : hd[0x0018, 0x0088].value.real,
            'AcquisitionDuration' : hd[0x0018, 0x9073].value.real,
            'RepetitionTime' : hd[0x0018,0x0080].value.real/1000,
#            'EchoTime' : hd[0x0018,0x0081].value.real/1000,
            'FlipAngle' : hd[0x0018, 0x1314].value.real,
            'pixdim_1' : hd[0x0028, 0x0030].value[0].real,
            'pixdim_2' : hd[0x0028, 0x0030].value[1].real,
            'pixdim_3' : hd[0x0018, 0x0050].value.real,
            'dim_1' : hd[0x0028, 0x0010].value.real,
            'dim_2' : hd[0x0028, 0x0011].value.real,
            'dim_3' : hd[0x2001, 0x1018].value.real,
            'dim_4' : hd[0x0020, 0x0105].value.real,
            'ParallelReductionFactorInPlane' : hd[0x2005, 0x140f][0][0x0018, 0x9069].value.real,
            'ParallelAcquisitionTechnique' : hd[0x2005, 0x140f][0][0x0018,0x9078].value,
#            'Units': hd[0x0028,0x1054].value,
            'EchoTime1':d_mag['EchoTime'],
            'EchoTime2':d_mag['EchoTime']+0.001, # EchoTime diff is 1ms (Peter)
            'IntendedFor': 'func/'+sub_folder+'_task-rest_bold.nii.gz'
            }   
        
        ## create the JSON file and save it in the func folder
        with open(dst_path_fmap+sub_folder+'_phasediff.json', 'w') as fp:
            json.dump(d_fp, fp, sort_keys = True, indent = 4,
                   ensure_ascii = False)   
        
        # save in dataframe
        df_dcm_fp = pd.DataFrame.from_dict(data=d_fp, orient='index').T
        df_dcm_fp.index = [df_dcm_fp.PatientsName.values[0].replace('-','_')]
        df_dcm_fp_all = pd.concat([df_dcm_fp_all, df_dcm_fp], sort=False)   

        

        
# sort dicom all
df_dcm_rsfmri_all = df_dcm_rsfmri_all.sort_index() 
df_dcm_rsfmri_all.insert(0, value=range(0,df_dcm_rsfmri_all.shape[0]), column='num_subjects')
df_dcm_t1_all = df_dcm_t1_all.sort_index() 
df_dcm_t1_all.insert(0, value=range(0,df_dcm_t1_all.shape[0]), column='num_subjects')
df_dcm_flair_all = df_dcm_flair_all.sort_index() 
df_dcm_flair_all.insert(0, value=range(0,df_dcm_flair_all.shape[0]), column='num_subjects')
df_dcm_fp_all = df_dcm_fp_all.sort_index() 
df_dcm_fp_all.insert(0, value=range(0,df_dcm_fp_all.shape[0]), column='num_subjects')
df_dcm_mag_all = df_dcm_mag_all.sort_index() 
df_dcm_mag_all.insert(0, value=range(0,df_dcm_mag_all.shape[0]), column='num_subjects')

# save
df_dcm_rsfmri_all.to_csv(img_path+'7T_dicom_rsfmri.csv')  
df_dcm_t1_all.to_csv(img_path+'7T_dicom_t1.csv')  
df_dcm_flair_all.to_csv(img_path+'7T_dicom_flair.csv')   
df_dcm_fp_all.to_csv(img_path+'7T_dicom_fieldmap.csv')  
df_dcm_mag_all.to_csv(img_path+'7T_dicom_magnitude.csv')  


print('>>> T1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
for item in df_dcm_t1_all:
    print(df_dcm_t1_all.groupby(item).count())

print('>>> FLAIR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
for item in df_dcm_flair_all:
    print(df_dcm_flair_all.groupby(item).count())
    
print('>>> rsfMRI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
for item in df_dcm_rsfmri_all:
    if item!='SliceTiming':
        print(df_dcm_rsfmri_all.groupby(item).count())

print('>>> fieldmap <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
for item in df_dcm_fp_all:
    print(df_dcm_fp_all.groupby(item).count())

print('>>> magnitude <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
for item in df_dcm_mag_all:
    print(df_dcm_mag_all.groupby(item).count())








print('>>>>>>>>>>>>>> The END')