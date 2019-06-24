############################################################
##
## useful parameters for seeded watershed
##
## to be used with 2-mars.py
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
##
## begin: first time point
## delta: delta/time interval between two time points (if one does not want
##        to deal with every single time point) (default = 1)
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
## Mars segmentation results will be stored in 'PATH_EMBRYO'/SEG/SEG_'EXP_SEG'
## (Mars is a seeded watershed procedure)
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
## mars parameters
##
######################################################################

# mars_begin = -1
# mars_end = -1

## ##### explanation #####
##
## interval of time points to be segmented by MARS method
## default is that only the time point 'begin' (see above)
## is segmented. More time points can be segmented only
## 'mars_begin' and 'mars_end' variables. Please note that
## a 'delta' (see above) step will take place between two processed time points.



# mars_intensity_transformation = 'Identity'
# mars_intensity_enhancement = None

## ##### explanation #####
##
## MARS method (nothing but a seeded watershed) may be applied on
## a transformed input image or on the original image (eg the result of
## the fusion step). This transformed imaged is made of a combination
## (by the maximum operator) of the transformed intensity image and a
## membrane-enhanced image.
## The transformed intensity image can be
## - None
## - Identity (the input image)
## - Normalization_to_u8 (a normalized version of the input image on 1 byte)
## The membrane-enhanced image can be
## - None
## - GACE: 'Global Automated Cell Extractor'
##
## This choice can be tuned with the two parameters
## - mars_intensity_transformation
## - mars_intensity_enhancement
##
## mars_intensity_transformation can be chosen in [None, 'Identity', 'Normalization_to_u8']
## mars_intensity_enhancement can be chosen in [None, 'GACE']
##



# mars_keep_reconstruction = True

## ##### explanation #####
##
## Previous tuning enables to perform the watershed segmentation on an image that is not
## the original image. If one does not want to keep such images, please turn the next variable to False



# mars_sigma_membrane = 0.9
# mars_hard_threshold = 1.0
# mars_sensitivity = 0.99
# mars_manual = False
# mars_manual_sigma = 15
# mars_sigma_TV = 3.6
# mars_sample = 0.2

## ##### explanation #####
##
## GACE parameters
## GACE method is made of three steps
## 1. membrane detection
## 2. membrane segmentation
## 3. tensor voting
##
## - membrane detection:
## this is the gaussian sigma that is used to compute image derivatives
## (in real units, a priori 0.9 um is a good choice for data like Patrick/Ralph/Aquila)
## mars_sigma_membrane=0.9
##
## - membrane segmentation:
## segmentation of the computed extrema is done by thresholding.
## Either mars_hard_thresholding is set to True, and a global hard threshold (mars_hard_threshold)
## is applied to the whole extrema image. This is not advised, and should be used only when the
## adaptative thresholding has failed.
## Or mars_hard_thresholding is set to False, meaning that an adaptative thresholding is computed
## by histogram fitting. This adaptative thresholding is governed by three parameters
## - mars_sensitivity:
##   membrane binarization parameter, /!\ if failure, one should enter in "manual" mode of the function
##   anisotropicHist via activation of 'manual' option
## - mars_manual:
##   By default, this parameter is set to False. If failure, (meaning that thresholds are very bad,
##   meaning that the binarized image is very bad), set this parameter to True and relaunch the
##   computation on the test image. If the method fails again, "play" with the value of manual_sigma...
##   and good luck.
## - mars_manual_sigma:
##   Axial histograms fitting initialization parameter for the computation of membrane image binarization
##   axial thresholds (this parameter is used if manual = True). One may need to test different values of
##   manual_sigma. We suggest to test values between 5 and 25 in case of initial failure. Good luck.
##
## - tensor voting:
## tensor voting is governed by 3 parameters
## - mars_sigma_TV:
##   parameter which defines the voting scale for membrane structures propagation by tensor voting method (real
##   coordinates). This parameter shoud be set between 3 um (little cells) and 4.5 um
##   (big gaps in the binarized membrane image)
## - mars_sample:
##   Parameter for tensor voting computation speed optimisation (do not touch if not bewared)
## - mars_sigma_LF:
##   Additional smoothing parameter for reconstructed image (in real coordinates).
##   It seems that the default value = 0.9 um is ok for standard use.
##



# watershed_seed_sigma = 0.6
# watershed_membrane_sigma = 0.15
# watershed_seed_hmin = 4

## ##### explanation #####
##
## watershed parameters
##
## the watershed segmentation has several steps
## 1. smoothing of initial image for seed extraction
##    mars_sigma1 (for back-compatibility)
##    watershed_seed_sigma
##    default value is 0.6 (real coordinates, ie 0.6 um)
## 2. smoothing of reconstructed image for image regularization prior to segmentation
##    mars_sigma2 (for back-compatibility)
##    watershed_membrane_sigma
##    default value is 0.15 (real coordinates, ie 0.15 um)
## 3. regional minima extraction
##    mars_h_min (for back-compatibility)
##    watershed_seed_hmin
##    default value is 4
##    Please note that this value may depend on the transformation applied on the input image, ie
##    h_min value should not be the same for the original image, or the image normalized into 1 byte.
##
