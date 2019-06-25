############################################################
##
## useful parameters for segmentation propagation
##
## to be used with 4-astec.py
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
delta = 1
raw_delay = 0

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



# EXP_FUSE = ''

## ##### explanation #####
##
## fusion directory name
##
## This fusion directory contains the images to be segmented
## default value is
## EXP_FUSE = 'RELEASE'
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



######################################################################
##
## astec parameters (new version, coming soon)
##
######################################################################



# astec_intensity_transformation = 'Identity'
# astec_intensity_enhancement = None
# astec_cell_normalization_min_method = 'cellinterior'
# astec_cell_normalization_max_method = 'cellborder'

## ##### explanation #####
##
## ASTEC method (nothing but a cell-based seeded watershed) may be applied on
## a transformed input image or on the original image (eg the result of
## the fusion step). This transformed imaged is made of a combination
## (by the maximum operator) of the transformed intensity image and a
## membrane-enhanced image.
## The transformed intensity image can be
## - None
## - Identity (the input image)
## - Normalization_to_u8 (a normalized version of the input image on 1 byte)
## - Cell_Normalization_to_u8 (a cell-based normalized version of the input image on 1 byte)
## The membrane-enhanced image can be
## - None
## - GACE: 'Global Automated Cell Extractor'
## - GLACE:
##
## This choice can be tuned with the two parameters
## - astec_intensity_transformation
## - astec_intensity_enhancement
##
## astec_intensity_transformation can be chosen in [None, 'Identity', 'Normalization_to_u8', 'Cell_Normalization_to_u8']
## astec_intensity_enhancement can be chosen in [None, 'GACE', 'GLACE']
##
## The 'Cell_Normalization_to_u8' intensity transformation method can be tuned with
## - astec_cell_normalization_min_method
## - astec_cell_normalization_min_method
## both variable are to be chosen in ['global', 'cell', 'cellborder', 'cellinterior', 'voxel']
## Choosing both of them as 'global' comes to choose 'Normalization_to_u8'
##



# astec_keep_reconstruction = True

## ##### explanation #####
##
## Previous tuning enables to perform the watershed segmentation on an image that is not
## the original image. If one does not want to keep such images, please turn the next variable to False



##########################################
##
## SEGMENTATION PROPAGATION PARAMETERS ###
##
##########################################



## Modules choice

astec_membrane_reconstruction_method = 0

## Membrane reconstruction module choice
## 0 for 'Classic' method
## 1 for 'Glace' method
## 2 for 'Gace' method
## If not set or set to 0, the input fused image is not processed for 
## membrane structures enhancement.
## If set to 1, the GLACE reconstruction method is going to be called
## If set to 2, the GACE reconstruction method is going to be called

astec_fusion_u8_method = 0

## Selection of the method which converts the fused image into a 8 bits
## images for the segmentation propagation. 
## If set to 0 (default), calling the historical "to_u8" method
## If set to 1, calling the mc-adhocFuse program which enhances the fused 
## image while converting it to u8 knowing the semgnetation propagation from
## previous time point

astec_flag_hybridation = False

## If set to True and if the membrane_reconstruction_method parameter is 
## provided and not equal to 0, then the reconstructed gray level image
## used for the segmentation propagation framework is goind to be an hybridation 
## between the original fused image and the result of image reconstruction by
## the specified method.

astec_keep_reconstruct_files = False 
## Set it to True in order to keep a copy
## Flag enabling to keep a copy of graylevel files provided to the watershed



## General parameters for segmentation propagation

astec_sigma1 = 0.6  		
## sigma 1 (0.6um) in real coordinates

astec_sigma2 = 0.15 		
## sigma 2 (0.15um) in real coordinates

astec_h_min_min = 4
## H min initialisation to ease correction

astec_h_min_max = 18   		
## H min initialisation to ease correction



## Glace Parameters (if astec_membrane_reconstruction_method is set to 1 or 2):
## membrane_renforcement

astec_sigma_membrane = 0.9
## membrane enhancement parameter (in real units, a
## priori 0.9 um is a good choice for data like 
## Patrick/Ralph/Aquila)

## anisotropicHist /!\ critical step
astec_sensitivity = 0.99  
## membrane binarization parameter, /!\ if failure,
## one should enter in "manual" mode of the function
## anisotropicHist via activation of 'manual' option

astec_manual = False     	
## By default, this parameter is set to False. If 
## failure, (meaning that thresholds are very bad, 
## meaning that the binarized image is very bad),
## set this parameter to True and relaunch the 
## computation on the test image. If the method fails
## again, "play" with the value of manual_sigma... 
## and good luck.

astec_manual_sigma = 15   
## Axial histograms fitting initialization parameter 
## for the computation of membrane image binarization
## axial thresholds (this parameter is used iif 
## manual = True).
## One may need to test different values of 
## manual_sigma. We suggest to test values between 5 and
## 25 in case of initial failure. Good luck.

astec_hard_thresholding = False 
## If the previous membrane threshold method 
## failed, one can force the thresholding with a
## "hard" threshold applied on the whole image. 
## To do so, this option must be set to True.

astec_hard_threshold = 1.0      
## If hard_thresholding = True, the enhanced 
## membranes image is thresholded using this 
## parameter (value 1 seems to be ok for 
## time-point t001 of Aquila embryo for example).



## Tensor voting framework

astec_sigma_TV = 3.6    
## parameter which defines the voting scale for membrane
## structures propagation by tensor voting method (real
## coordinates). 
## This parameter shoud be set between 3 um (little cells)
## and 4.5 um(big gaps in the binarized membrane image)

astec_sigma_LF = 0.9    
## Smoothing parameter for reconstructed image (in real
## coordinates). It seems that the default value = 0.9 um
## is ok for classic use.

astec_sample = 0.2      
## Parameter for tensor voting computation speed 
## optimisation (do not touch if not bewared)

astec_rayon_dil = 3.6   
## dilatation ray for propagated ROI from time t to t+1
## (default: 3.6, in real coordinates) 



## Fused image conversion parameters (if astec_fusion_u8_method is set to 1)

astec_min_percentile = 0.01   
## mc-adhocFuse parameter of type %f (default: 0.01)

astec_max_percentile = 0.99   
## mc-adhocFuse parameter of type %f (default: 0.99)

astec_min_method = 'cellinterior'
## mc-adhocFuse param. (default: 'cellinterior')
## taken in global|cell|cellborder|cellinterior|voxel

astec_max_method = 'cellborder'
## mc-adhocFuse parameter (default: 'cellborder')
## taken in global|cell|cellborder|cellinterior|voxel

astec_sigma_hybridation = 5.0 
## mc-adhocFuse parameter of type %f (default: 5.0)



## Default parameters (for classical use, default values should not be changed)

astec_RadiusOpening = 20 		
## (using the approximation of a sphere of radius 20
## voxels as a structuring element)

astec_Thau = 25 				
## s(c)=h2+(c).N2(c) >t identical

astec_MinVolume = 1000 		
## Miss Suppressing cells with to small volumes (not
## yet in supdata)

astec_VolumeRatioBigger = 0.5 
## If a cell in St+1 is at least 50% bigger than its
## progeny in St+1,

astec_VolumeRatioSmaller = 0.1
## Cells in St+1 that are 10% or smaller than their
## equivalent in St+1 are tagged for correction

astec_MorphosnakeIterations = 10 
## Then, an active contour algorithm is applied
## using the dilated shape of c, obtained by 
## iterating 10 times

astec_NIterations = 200 		
## The algorithm is applied up to stability 
## (at th voxels) or after n iterations 
## (in our case th = 103 and n = 200).

astec_DeltaVoxels = 10**3  	
## y (at th voxels)

astec_Volum_Min_No_Seed = 100 
## Then, if the volume of c is greater than 100 
## voxels (2.7 um3)

astec_nb_proc = 10 			
## Number of processor ...

astec_nb_proc_ace = 7   		
## number of processors for ACE (7 is recommanded)



