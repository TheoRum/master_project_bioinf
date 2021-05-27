%{
 | Title:
     plotting the connectograms
     
 | Date:
     2021-04-10
     
 | Author(s):
     Theodor Rumetshofer

 | Description:    
     This script specifies the input data for the GIFT function

 | List of functions:
     No user defined functions are used in the program.

 | List of "non standard" modules:
     None

 | Procedure:
     For CPAC, fMRIprep and  the differences of the correlation matrices 
     were conectograms generated using the GIFT function icatb_plot_connectogram.m

%}


%% session 1 (CPAC)

load('/Volumes/ID1036/rsfMRI_7T/gift210430/output_5/network_summary_top_s1/connectogram_s1_cpac.mat')

% thresh the conn_threshold=1.5
icatb_plot_connectogram([], comp_network_names, 'C',C, 'image_file_names',image_file_names, 'convert_to_zscores', 'yes', 'threshold', 1.96, 'title', '', 'image_values', 'positive', 'template_file', template_file, 'CLIM', [-0.4, 0.6], 'conn_threshold', 0.15, 'imwidth', 0.15, 'colorbar_label', 'Correlation')
print -r300 -dpng /Users/theo/Dropbox/work/LundUniversity/publications_abstracts/SLE_rsfMRI/Master_report/figures/connectogram_connthresh0.15_cpac.png
close all

% no conn_threshold
icatb_plot_connectogram([], comp_network_names, 'C',C, 'image_file_names',image_file_names, 'convert_to_zscores', 'yes', 'threshold', 1.96, 'title', '', 'image_values', 'positive', 'template_file', template_file, 'CLIM', [-0.4, 0.6], 'imwidth', 0.15, 'colorbar_label', 'Correlation')
print -r300 -dpng /Users/theo/Dropbox/work/LundUniversity/publications_abstracts/SLE_rsfMRI/Master_report/figures/connectogram_cpac.png
close all

clear all

%% session 2 (fMRIprep)

load('/Volumes/ID1036/rsfMRI_7T/gift210430/output_5/network_summary_top_s2/connectogram_s2_fmriprep.mat')

% thresh the conn_threshold=1.5
icatb_plot_connectogram([], comp_network_names, 'C',C, 'image_file_names',image_file_names, 'convert_to_zscores', 'yes', 'threshold', 1.96, 'title', '', 'image_values', 'positive', 'template_file', template_file, 'CLIM', [-0.4, 0.6], 'conn_threshold', 0.15, 'imwidth', 0.15, 'colorbar_label', 'Correlation')
print -r300 -dpng /Users/theo/Dropbox/work/LundUniversity/publications_abstracts/SLE_rsfMRI/Master_report/figures/connectogram_connthresh0.15_fmriprep.png
close all

% no conn_threshold
icatb_plot_connectogram([], comp_network_names, 'C',C, 'image_file_names',image_file_names, 'convert_to_zscores', 'yes', 'threshold', 1.96, 'title', '', 'image_values', 'positive', 'template_file', template_file, 'CLIM', [-0.4, 0.6], 'imwidth', 0.15, 'colorbar_label', 'Correlation')
print -r300 -dpng /Users/theo/Dropbox/work/LundUniversity/publications_abstracts/SLE_rsfMRI/Master_report/figures/connectogram_fmriprep.png
close all

clear all

%% session diff (fMRIprep-CPAC)

load('/Volumes/ID1036/rsfMRI_7T/gift210430/output_5/network_summary_top_s_diff/connectogram_s_diff.mat')

% thresh the conn_threshold=1.5
icatb_plot_connectogram([], comp_network_names, 'C',C, 'image_file_names',image_file_names, 'convert_to_zscores', 'yes', 'threshold', 1.96, 'title', '', 'image_values', 'positive', 'template_file', template_file, 'CLIM', [-0.4, 0.5], 'conn_threshold', 0.15, 'imwidth', 0.15, 'colorbar_label', 'Correlation')
print -r300 -dpng /Users/theo/Dropbox/work/LundUniversity/publications_abstracts/SLE_rsfMRI/Master_report/figures/connectogram_connthresh0.15_diff_f-c.png
close all

% thresh the conn_threshold=1.5
icatb_plot_connectogram([], comp_network_names, 'C',C, 'image_file_names',image_file_names, 'convert_to_zscores', 'yes', 'threshold', 1.96, 'title', '', 'image_values', 'positive', 'template_file', template_file, 'CLIM', [-0.4, 0.5], 'imwidth', 0.15, 'colorbar_label', 'Correlation')
print -r300 -dpng /Users/theo/Dropbox/work/LundUniversity/publications_abstracts/SLE_rsfMRI/Master_report/figures/connectogram_diff_f-c.png
close all

