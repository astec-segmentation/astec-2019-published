############################################################
##
## useful parameters for segmentation correction
##
## to be used with 3-manualcorrection-parameters.py
##
############################################################


# PATH_EMBRYO = ''

## ##### explanation #####
##
## path to the embryo data
## e.g. '/media/DATA/171107-Karine-St8'
## if not present, the current directory is used
##



# EN = ''

## ##### explanation #####
##
## Embryo Name
## CRBM naming format is YYMMDD-SaintOfTheDays-Stage
## eg: '171107-Karine-St8'
## (automatically extracted from PATH_EMBRYO if not provided)
##



begin = 0
# delta = 1

## ##### explanation #####
##
## time point to be processed
## to process several time points, use 'mars_begin' and 'mars_end'
## warning: the same correction file parameter will be applied to all the time points 
##
## begin: first time point
## delta: delta/time interval between two time points (if one does not want
##        to deal with every single time point) (default = 1)
##



# EXP_SEG = ''

## ##### explanation #####
##
## Path to both the image to be corrected and the result image to be saved
##
## Segmentation (Astec) results will be stored in 'PATH_EMBRYO'/SEG/SEG_'EXP_SEG'
## default value is
## EXP_SEG = 'RELEASE'
##



# result_image_suffix = 'mha'
# default_image_suffix = 'mha'

## ##### explanation #####
##
## defines the output image format
##
## result_image_suffix is for all output images
## default_image_suffix is for all images (including auxiliary ones)
## default are 'inr'
##
## 'mha' should be prefered (readable by fiji)
##



######################################################################
##
## manual correction parameters
##
######################################################################



# mars_begin = -1
# mars_end = -1

## ##### explanation #####
##
## interval of time points to be corrected
## default is that only the time point 'begin' (see above)
## is segmented. More time points can be corrected only
## 'mars_begin' and 'mars_end' variables. Please note that
## a 'delta' (see above) step will take place between two processed time points.



# mancor_input_seg_file = ''
# mancor_output_seg_file = ''

## ##### explanation #####
##
## defines the input/ouput file names (to be used when correcting
## other files than the 2-mars.py output file)
## 
## mancor_input_seg_file is the segmentation file to be corrected
##   if not provided, then looking for the output of MARS
## mancor_output_seg_file is the corrected segmentation file                    



mancor_mapping_file = ''

## ##### explanation #####
##
## path to mapping file for manual correction of a 
## segmentation (ie label) image. See above the syntax of this file.
## - 1 line per label association
## - background label has value 1
## - the character '#' denotes commented lines 

'''
# EXAMPLE OF mancor_mapping_file CONTENT:
# here the input label 8 will be mapped with new value 7, etc...
8 7
9 2  
4 64 
29 23
# ... etc ...
# background labels
30 1 
89 1 
'''

