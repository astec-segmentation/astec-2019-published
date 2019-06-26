############################################################
##
## useful parameters for segmentation post-correction
##
## to be used with 5-postcorrection.py 
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
end = 0
# delta = 1
# raw_delay = 0

## ##### explanation #####
##
## time points to be processed
##
## begin: first time point
## end: last time point
## delta: delta/time interval between two time points (if one does not want
##        to deal with every single time point) (default = 1)
##
## For the fusion step, not giving the 'begin' and 'end' values
## causes the fusion of all found data in raw data directories
## according these 4 directories are different
##
## When testing or tuning parameters, it is advised to set 'begin' and 'end'
## at the same value 
##



# EXP_SEG = ''

## ##### explanation #####
##
## Segmentation results will be stored in 'PATH_EMBRYO'/SEG/SEG_'EXP_SEG'
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



#################################
##
## POST-CORRECTION PARAMETERS ###
##
#################################

postcor_Volume_Threshold=10000 	
                      # volume low threshold for final cells 
					  # (in voxels) 
postcor_Soon=True 				
                      # True if the cell life span has to be 
					  # taken into account
postcor_ShortLifespan=25 		
                      # (length < SL time points, in our case 
					  # SL = 10)
postcor_PearsonThreshold=0.9; 	
                      # If they are anticorrelated (Pearson 
					  # correlation under -0.9)



