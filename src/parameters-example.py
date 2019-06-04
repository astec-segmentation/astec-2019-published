

# File for parameters and nomenclature settings

# ##### explanation #####
#
# path to the embryo data
# e.g. '/media/DATA/171107-Karine-St8'
# Can also be set by providing option '-e' in the command line
# if not present, the current directory is used
#

# PATH_EMBRYO = ''


# ##### explanation #####
#
# Embryo Name
# CRBM naming format is YYMMDD-SaintOfTheDays-Stage
# eg: '171107-Karine-St8'
# (automatically extracted from PATH_EMBRYO if not provided)
#

# EN = ''


######################################################################
#
# paths
#
######################################################################

# ##### explanation #####
#
# allows to specify paths to the raw (eq acquisition) data
# the 4 paths are built as follows
# 'PATH_EMBRYO'/'DIR_RAWDATA'/'DIR_LEFTCAM_STACKZERO'
# 'PATH_EMBRYO'/'DIR_RAWDATA'/'DIR_RIGHTCAM_STACKZERO'
# 'PATH_EMBRYO'/'DIR_RAWDATA'/'DIR_LEFTCAM_STACKONE'
# 'PATH_EMBRYO'/'DIR_RAWDATA/'DIR_RIGHTCAM_STACKONE'
#
# default values are
#
# DIR_RAWDATA = 'RAWDATA'
# DIR_LEFTCAM_STACKZERO = 'LC/Stack0000'
# DIR_RIGHTCAM_STACKZERO = 'RC/Stack0000'
# DIR_LEFTCAM_STACKONE = 'LC/Stack0001'
# DIR_RIGHTCAM_STACKONE = 'RC/Stack0001'
#
# for multi-channel fusion, paths to the raw data have also to be specified
# with the variables (X = 2 or 3) and the same path construction
# eg: 'PATH_EMBRYO'/'DIR_RAWDATA_CHANNEL_X'/'DIR_LEFTCAM_STACKZERO_CHANNEL_X'
# DIR_RAWDATA_CHANNEL_X
# DIR_LEFTCAM_STACKZERO_CHANNEL_X
# DIR_RIGHTCAM_STACKZERO_CHANNEL_X
# DIR_LEFTCAM_STACKONE_CHANNEL_X
# DIR_RIGHTCAM_STACKONE_CHANNEL_X
#
# if DIR_RAWDATA_CHANNEL_X is not given, it is replaced by DIR_RAWDATA or its default value
# if any of the four path DIR_dirCAM_STACKx_CHANNEL_X is not given, it is also replaced
# by the main channel variable (DIR_dirCAM_STACKx) or its default value
#

# DIR_RAWDATA = ''
# DIR_LEFTCAM_STACKZERO = ''
# DIR_RIGHTCAM_STACKZERO = ''
# DIR_LEFTCAM_STACKONE = ''
# DIR_RIGHTCAM_STACKONE = ''

# DIR_RAWDATA_CHANNEL_2 = ''
# DIR_LEFTCAM_STACKZERO_CHANNEL_2 = ''
# DIR_RIGHTCAM_STACKZERO_CHANNEL_2 = ''
# DIR_LEFTCAM_STACKONE_CHANNEL_2 = ''
# DIR_RIGHTCAM_STACKONE_CHANNEL_2 = ''

# DIR_RAWDATA_CHANNEL_3 = ''
# DIR_LEFTCAM_STACKZERO_CHANNEL_3 = ''
# DIR_RIGHTCAM_STACKZERO_CHANNEL_3 = ''
# DIR_LEFTCAM_STACKONE_CHANNEL_3 = ''
# DIR_RIGHTCAM_STACKONE_CHANNEL_3 = ''

# ##### explanation #####
#
# allows to specify prefix name of the raw data
# may be useful when several camera acquisitions are stored into
# the same directory. Eg, if one has set
#   DIR_RAWDATA='raw_data'
#   DIR_LEFTCAM_STACKZERO = 'Stack_0_Channel_0'
#   DIR_RIGHTCAM_STACKZERO = 'Stack_0_Channel_0'
#   DIR_LEFTCAM_STACKONE = 'Stack_1_Channel_0'
#   DIR_RIGHTCAM_STACKONE = 'Stack_1_Channel_0'
# meaning that stacks #0 of both left and right cameras are stored into the
# same directory 'raw_data/Stack_0_Channel_0'
#
# File names for the left camera acquisition are then built with
#   'PATH_EMBRYO'/'DIR_RAWDATA'/'DIR_TO_LEFTCAM_STACKZERO'/'acquisition_leftcam_image_prefix'+'TIME'
#   'PATH_EMBRYO'/'DIR_RAWDATA'/'DIR_TO_LEFTCAM_STACKONE'/'acquisition_leftcam_image_prefix'+'TIME'
# where 'TIME' is a 3-digit number (left-filled with 0's). The same for the right camera
# acquisition files. Default values are
#
# acquisition_leftcam_image_prefix = 'Cam_Left_00'
# acquisition_rightcam_image_prefix = 'Cam_Right_00'
#

# acquisition_leftcam_image_prefix = ''
# acquisition_rightcam_image_prefix = ''


# ##### explanation #####
#
# allows to specify path to the results
#
# Workspace for the considered step; this workspace
# is included in the repository corresponding to the
# ASTEC step ([FUSE|SEG|POST]) and prefixed by this
# repository name (eg. FUSE_EXP_NO_CROP or SEG_RELEASE)
# ----> use replaceDIR_EXP_FUSE(path_fuse_exp[...], EXP_FUSE)
# method from <astec-package>/nomenclature.py to get the right
# path for the given fusion experience
#


# ##### explanation #####
#
# Fusion results will be stored in 'PATH_EMBRYO'/FUSE/FUSE_'EXP_FUSE'
# default value is
# EXP_FUSE = 'RELEASE'
#
# for multi-channel fusion, a fusion result per additional channel has to be given
# EXP_FUSE_CHANNEL_2 = ''
# EXP_FUSE_CHANNEL_3 = ''
#

# EXP_FUSE = ''
# EXP_FUSE_CHANNEL_2 = ''
# EXP_FUSE_CHANNEL_3 = ''

#
#
#


# EXP_INTRAREG = ''


# ##### explanation #####
#
# Mars segmentation results will be stored in 'PATH_EMBRYO'/SEG/SEG_'EXP_MARS'
# default value is
# EXP_MARS = 'RELEASE'
#

# EXP_MARS=''


# ##### explanation #####
#
# Segmentation (Astec) results will be stored in 'PATH_EMBRYO'/SEG/SEG_'EXP_SEG'
# default value is
# EXP_SEG = 'RELEASE'
#

# EXP_SEG = ''


# ##### explanation #####
#
# Post-segmentation results will be stored in 'PATH_EMBRYO'/POST/POST_'EXP_POST'
# default value is
# EXP_POST = 'RELEASE'
#

# EXP_POST=''


#
#
#

# TIME=''		# Time-point of an embryo snapshot
# TIMEREF=''		# For registration, time-point of reference snapshot
# TIMEFLO=''		# For registration, time-point of floating snapshot


######################################################################
#
# general control parameters
#
######################################################################

# ##### explanation #####
#
# defines the acquisition time points to be processed
# begin: first time point
# end: last time point
# delta: Delta between two time points (if one does not want
#        to deal with every single time point) (default = 1)
#
# For the fusion step, not giving the 'begin' and 'end' values
# causes the fusion of all found data in raw data directories
# according these 4 directories are different
#

# begin = 0

# end = 0

# delta = 1


# ##### explanation #####
#
# defines the output image format
#
# RESULT_IMAGE_SUFFIX_FUSE is for the fusion result image(s)
# result_image_suffix is for all output images
# default_image_suffix is for all images (including auxiliary ones)
# default are 'inr'
#

# RESULT_IMAGE_SUFFIX_FUSE = 'inr'
# result_image_suffix = 'inr'
# default_image_suffix = 'inr'

######################################################################
#
# raw data definition
#
######################################################################

# ##### explanation #####
#
# raw_ori: image orientation
#   if im2 angle - im1 angle < 0 => raw_ori = 'right'
# raw_mirrors: depends on the acquisition protocol, value
# 	can be set to True or False
# 	- the standard value of this parameter is False
# 	- in case of axial symmetry between left
# 	  and right cameras, then set to True
# raw_resolution: acquisition voxel size
#   e.g. raw_resolution = (.21, .21, 1.)
# raw_delay: increment to to be added to the time values
#   (values in range [begin,end]) when generating the
#   fused image names (default is 0)
#

raw_ori = 'left'

raw_mirrors = False

raw_resolution = (.17, .17, 1.)

raw_delay = 0


######################################################################
#
# fusion parameters
#
######################################################################

# ##### explanation #####
#
# the fusion of the 4 acquisitions follows a number of steps
#
# 1. Optionally, a slit line correction.
#    Some Y lines may appear brighter in the acquisition and causes artifacts in the
#    reconstructed (ie fused) image.
# 2. a change of resolution in the X and Y directions only (Z remains unchanged)
#    it allows to decrease the data volume if the new pixel size is larger
#    than the acquisition one
# 3. Optionally, a crop of the resampled acquisitions
#    it allows to decrease the volume of data
#    the crop is based on the analysis of a MIP view (in the Z direction) of
#    the volume
# 4. Optionally, a mirroring of the 'right' image
#    depends of the 'raw_mirrors' value (see supra)
# 5. Linear registration of the 3 last images on the first one (considered as the reference)
#    The reference image is resampled again, to get an isotropic voxel
#    (same voxel size in the 3 directions: X, Y, Z)
# 6. Linear combination of images, weighted by an ad-hoc function
# 7. Crop of the fused image
#    still based on the analysis of a MIP view (in the Z direction)
#

#
# step 1. parameters (also used in step 5.)
# Isotropic resolution of the final fused image
#

# acquisition_slit_line_correction = True

#
# step 2. parameters (also used in step 5.)
# Isotropic resolution of the final fused image
#

target_resolution = .3


#
# step 3. parameters
#
# raw_crop: if False, then the acquisition images are not cropped
#   (default is True)
# raw_margin_x_0: parameter for margin of the bounding box computed
#   for the cropping of the  raw acquisition image in 'left' X direction
#   (default is 40)
# raw_margin_x_1: parameter for margin of the bounding box computed
#   for the cropping of the  raw acquisition image in 'right' X direction
#   (default is 40)
# raw_margin_y_0: parameter for margin of the bounding box computed
#   for the cropping of the  raw acquisition image in 'left' Y direction
#   (default is 40)
# raw_margin_y_1: parameter for margin of the bounding box computed
#   for the cropping of the  raw acquisition image in 'right' Y direction
#   (default is 40)
#

# raw_crop = True
# raw_margin_x_0 = 40
# raw_margin_x_1 = 40
# raw_margin_y_0 = 40
# raw_margin_y_1 = 40

#
# step 5. parameters
#
# fusion_registration_pyramid_highest_level: highest level of the pyramid image for registration
#   registration is done hierarchically with a pyramid of images at each level, image dimensions are divided by 2.
#   'fusion_registration_pyramid_highest_level = 6' means that registration starts with images whose dimensions
#   are 1/64th of the original image
# fusion_registration_pyramid_lowest_level: lowest level of the pyramid image for registration
#

# fusion_preregistration_compute_registration = False
# fusion_preregistration_transformation_type = 'translation'
# fusion_preregistration_transformation_estimation_type = 'wlts'
# fusion_preregistration_lts_fraction = 0.55
# fusion_preregistration_pyramid_highest_level = 6
# fusion_preregistration_pyramid_lowest_level = 3
# fusion_preregistration_normalization = True

# fusion_registration_compute_registration = True
# fusion_registration_transformation_type = 'affine'
# fusion_registration_transformation_estimation_type = 'wlts'
# fusion_registration_lts_fraction = 0.55
# fusion_registration_pyramid_highest_level = 6
# fusion_registration_pyramid_lowest_level = 3
# fusion_registration_normalization = True

#
# step 7. parameters
#
# fusion_crop: if False, then the acquisition images are not cropped
#   (default is True)
# fusion_margin_x_0: parameter for margin of the bounding box computed
#   for the cropping of the  raw acquisition image in 'left' X direction
#   (default is 40)
# fusion_margin_x_1: parameter for margin of the bounding box computed
#   for the cropping of the  raw acquisition image in 'right' X direction
#   (default is 40)
# fusion_margin_y_0: parameter for margin of the bounding box computed
#   for the cropping of the  raw acquisition image in 'left' Y direction
#   (default is 40)
# fusion_margin_y_1: parameter for margin of the bounding box computed
#   for the cropping of the  raw acquisition image in 'right' Y direction
#   (default is 40)
#

# fusion_crop = True
# fusion_margin_x_0 = 40
# fusion_margin_x_1 = 40
# fusion_margin_y_0 = 40
# fusion_margin_y_1 = 40


######################################################################
#
# intra-registration parameters
#
######################################################################

# ##### explanation #####
#
# Intra-registration results are stored into 'PATH_EMBRYO'/INTRAREG/INTRAREG_'EXP_INTRAREG'
#


#
# parameters for the co-registration of successive fused images
#

# intra_registration_compute_registration = True
# intra_registration_transformation_type = 'rigid'
# intra_registration_transformation_estimation_type = 'wlts'
# intra_registration_lts_fraction = 0.55
# intra_registration_pyramid_highest_level = 6
# intra_registration_pyramid_lowest_level = 3
# intra_registration_normalization = True

#
# parameters for template building
# the 'intra_registration_resolution' parameter gives the resulting (isotropic) voxel size
# (as the 'target_resolution' gives the voxel size of the fused images). However, for
# visualization purposes, it may be indicated to have a larger voxel size (hence the 0.6
# instead the 0.3)
#
# The template is built so that the useful information of all resampled images fits into it.
# Useful information can be issued from either the fused sequence, the segmentation sequence or
# the post-segmentation sequence. It can be indicated by the 'intra_registration_template_type'
# parameter:
#  - intra_registration_template_type = 'FUSION' | 'SEGMENTATION' | 'POST-SEGMENTATION'
# By default, the useful information is the whole image, which may yield a huge template image.
#
# Giving a threshold with the 'intra_registration_template_type', only points above the threshold
# are considered to be included in the template after resampling, this allows to reduce the template.
# According the background value is either 0 or 1 in both the segmentation and the post-segmentation
# sequences, setting this threshold to 2 for these sequences allows to keep the entire embryo in the
# resampled/reconstructed sequence.
#
# In addition, a margin can be given for a more comfortable visualization. By default, it is 0 when only
# fusion images are used, and 10 if either segmentation or post-segmentation images are also used
#

# intra_registration_reference_index = None
# intra_registration_template_type = "FUSION"
# intra_registration_template_threshold = None
# intra_registration_resolution = 0.6
# intra_registration_margin = None

#
# outputs
# 1. resampled images (in the 'template' geometry)
#    If required, it assumes that segmentation and post-segmentation images exists
#    and are respectively located in 'PATH_EMBRYO'/SEG/SEG_'EXP_SEG' and
#    'PATH_EMBRYO'/POST/POST_'EXP_POST'
# 2. movies (ie 2D+t images)
#    each 'intra_registration_[xy|xz|yz]_movie_fusion_images
#

# intra_registration_sigma_segmentation_images = 1.0

# intra_registration_resample_fusion_images = True
# intra_registration_resample_segmentation_images = False
# intra_registration_resample_post_segmentation_images = False

# intra_registration_movie_fusion_images = True
# intra_registration_movie_segmentation_images = False
# intra_registration_movie_post_segmentation_images = False

# intra_registration_xy_movie_fusion_images = [];
# intra_registration_xz_movie_fusion_images = [];
# intra_registration_yz_movie_fusion_images = [];

# intra_registration_xy_movie_segmentation_images = [];
# intra_registration_xz_movie_segmentation_images = [];
# intra_registration_yz_movie_segmentation_images = [];

# intra_registration_xy_movie_post_segmentation_images = [];
# intra_registration_xz_movie_post_segmentation_images = [];
# intra_registration_yz_movie_post_segmentation_images = [];

######################################################################
#
# mars parameters
#
######################################################################

# ##### explanation #####
#
# interval of time points to be segmented by MARS method
# default is that only the time point 'begin' (see above)
# is segmented. More time points can be segmented only
# 'mars_begin' and 'mars_end' variables. Please note that
# a 'delta' (see above) step will take place between two processed time points.

# mars_begin = -1
# mars_end = -1

# ##### explanation #####
#
# MARS method (nothing but a seeded watershed) may be applied on
# a transformed input image or on the original image (eg the result of
# the fusion step). This transformed imaged is made of a combination
# (by the maximum operator) of the transformed intensity image and a
# membrane-enhanced image.
# The transformed intensity image can be
# - None
# - Identity (the input image)
# - Normalization_to_u8 (a normalized version of the input image on 1 byte)
# The membrane-enhanced image can be
# - None
# - GACE: 'Global Automated Cell Extractor'
#
# This choice can be tuned with the two parameters
# - mars_intensity_transformation
# - mars_intensity_enhancement
#
# mars_intensity_transformation can be chosen in [None, 'Identity', 'Normalization_to_u8']
# mars_intensity_enhancement can be chosen in [None, 'GACE']
#

# mars_intensity_transformation = 'Identity'
# mars_intensity_enhancement = None


# ##### explanation #####
#
# Previous tuning enables to perform the watershed segmentation on an image that is not
# the original image. If one does not want to keep such images, please turn the next variable to False

# mars_keep_reconstruction = True

# ##### explanation #####
#
# GACE parameters
# GACE method is made of three steps
# 1. membrane detection
# 2. membrane segmentation
# 3. tensor voting
#
# - membrane detection:
# this is the gaussian sigma that is used to compute image derivatives
# (in real units, a priori 0.9 um is a good choice for data like Patrick/Ralph/Aquila)
# mars_sigma_membrane=0.9
#
# - membrane segmentation:
# segmentation of the computed extrema is done by thresholding.
# Either mars_hard_thresholding is set to True, and a global hard threshold (mars_hard_threshold)
# is applied to the whole extrema image. This is not advised, and should be used only when the
# adaptative thresholding has failed.
# Or mars_hard_thresholding is set to False, meaning that an adaptative thresholding is computed
# by histogram fitting. This adaptative thresholding is governed by three parameters
# - mars_sensitivity:
#   membrane binarization parameter, /!\ if failure, one should enter in "manual" mode of the function
#   anisotropicHist via activation of 'manual' option
# - mars_manual:
#   By default, this parameter is set to False. If failure, (meaning that thresholds are very bad,
#   meaning that the binarized image is very bad), set this parameter to True and relaunch the
#   computation on the test image. If the method fails again, "play" with the value of manual_sigma...
#   and good luck.
# - mars_manual_sigma:
#   Axial histograms fitting initialization parameter for the computation of membrane image binarization
#   axial thresholds (this parameter is used if manual = True). One may need to test different values of
#   manual_sigma. We suggest to test values between 5 and 25 in case of initial failure. Good luck.
#
# - tensor voting:
# tensor voting is governed by 3 parameters
# - mars_sigma_TV:
#   parameter which defines the voting scale for membrane structures propagation by tensor voting method (real
#   coordinates). This parameter shoud be set between 3 um (little cells) and 4.5 um
#   (big gaps in the binarized membrane image)
# - mars_sample:
#   Parameter for tensor voting computation speed optimisation (do not touch if not bewared)
# - mars_sigma_LF:
#   Additional smoothing parameter for reconstructed image (in real coordinates).
#   It seems that the default value = 0.9 um is ok for standard use.
#

# mars_sigma_membrane = 0.9
# mars_hard_threshold = 1.0
# mars_sensitivity = 0.99
# mars_manual = False
# mars_manual_sigma = 15
# mars_sigma_TV = 3.6
# mars_sample = 0.2


# ##### explanation #####
#
# watershed parameters
#
# the watershed segmentation has several steps
# 1. smoothing of initial image for seed extraction
#    mars_sigma1 (for back-compatibility)
#    watershed_seed_sigma
#    default value is 0.6 (real coordinates, ie 0.6 um)
# 2. smoothing of reconstructed image for image regularization prior to segmentation
#    mars_sigma2 (for back-compatibility)
#    watershed_membrane_sigma
#    default value is 0.15 (real coordinates, ie 0.15 um)
# 3. regional minima extraction
#    mars_h_min (for back-compatibility)
#    watershed_seed_hmin
#    default value is 4
#    Please note that this value may depend on the transformation applied on the input image, ie
#    h_min value should not be the same for the original image, or the image normalized into 1 byte.
#

# watershed_seed_sigma = 0.6
# watershed_membrane_sigma = 0.15
# watershed_seed_hmin = 4


######################################################################
#
#
#
######################################################################

####################################
### MANUAL CORRECTION PARAMETERS ###
####################################

mancor_input_seg_file=''
                       # segmentation file to be manually corrected
					   # If not provided, then looking for the output of MARS
mancor_output_seg_file=''
mancor_mapping_file='' 
                       # path to mapping file for manual correction of the 
					   # mars segmentation. See above the syntax of this file.
					   # - 1 line per label association
					   # - background label has value 1
					   # - the character '#' denotes commented lines 
					   # See file "mapping.txt" in the astec project to get an
					   # example
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


######################################################################
#
# astec parameters
#
######################################################################


# ##### explanation #####
#
# ASTEC method (nothing but a cell-based seeded watershed) may be applied on
# a transformed input image or on the original image (eg the result of
# the fusion step). This transformed imaged is made of a combination
# (by the maximum operator) of the transformed intensity image and a
# membrane-enhanced image.
# The transformed intensity image can be
# - None
# - Identity (the input image)
# - Normalization_to_u8 (a normalized version of the input image on 1 byte)
# - Cell_Normalization_to_u8 (a cell-based normalized version of the input image on 1 byte)
# The membrane-enhanced image can be
# - None
# - GACE: 'Global Automated Cell Extractor'
# - GLACE:
#
# This choice can be tuned with the two parameters
# - astec_intensity_transformation
# - astec_intensity_enhancement
#
# astec_intensity_transformation can be chosen in [None, 'Identity', 'Normalization_to_u8', 'Cell_Normalization_to_u8']
# astec_intensity_enhancement can be chosen in [None, 'GACE', 'GLACE']
#
# The 'Cell_Normalization_to_u8' intensity transformation method can be tuned with
# - astec_cell_normalization_min_method
# - astec_cell_normalization_min_method
# both variable are to be chosen in ['global', 'cell', 'cellborder', 'cellinterior', 'voxel']
# Choosing both of them as 'global' comes to choose 'Normalization_to_u8'

# astec_intensity_transformation = 'Identity'
# astec_intensity_enhancement = None
# astec_cell_normalization_min_method = 'cellinterior'
# astec_cell_normalization_max_method = 'cellborder'

# ##### explanation #####
#
# Previous tuning enables to perform the watershed segmentation on an image that is not
# the original image. If one does not want to keep such images, please turn the next variable to False

# astec_keep_reconstruction = True




###########################################
### SEGMENTATION PROPAGATION PARAMETERS ###
###########################################


# Modules choice

astec_membrane_reconstruction_method=0
# Membrane reconstruction module choice
# 0 for 'Classic' method
# 1 for 'Glace' method
# 2 for 'Gace' method
# If not set or set to 0, the input fused image is not processed for 
# membrane structures enhancement.
# If set to 1, the GLACE reconstruction method is going to be called
# If set to 2, the GACE reconstruction method is going to be called

astec_fusion_u8_method=0
# Selection of the method which converts the fused image into a 8 bits
# images for the segmentation propagation. 
# If set to 0 (default), calling the historical "to_u8" method
# If set to 1, calling the mc-adhocFuse program which enhances the fused 
# image while converting it to u8 knowing the semgnetation propagation from
# previous time point

astec_flag_hybridation=False
# If set to True and if the membrane_reconstruction_method parameter is 
# provided and not equal to 0, then the reconstructed gray level image
# used for the segmentation propagation framework is goind to be an hybridation 
# between the original fused image and the result of image reconstruction by
# the specified method.

astec_keep_reconstruct_files=False 
# Set it to True in order to keep a copy
# Flag enabling to keep a copy of graylevel files provided to the watershed



# General parameters for segmentation propagation
astec_sigma1 = 0.6  		
                            # sigma 1 (0.6um) in real coordinates
astec_sigma2 = 0.15 		
                            # sigma 2 (0.15um) in real coordinates
astec_h_min_min = 4
                            # H min initialisation to ease correction
astec_h_min_max = 18   		
                            # H min initialisation to ease correction

# Glace Parameters (if astec_membrane_reconstruction_method is set to 1 or 2):
# membrane_renforcement
astec_sigma_membrane=0.9
                        # membrane enhancement parameter (in real units, a
						# priori 0.9 um is a good choice for data like 
						# Patrick/Ralph/Aquila)
# anisotropicHist /!\ critical step
astec_sensitivity=0.99  
                        # membrane binarization parameter, /!\ if failure,
						# one should enter in "manual" mode of the function
						# anisotropicHist via activation of 'manual' option

astec_manual=False     	
                        # By default, this parameter is set to False. If 
						# failure, (meaning that thresholds are very bad, 
						# meaning that the binarized image is very bad),
				 		# set this parameter to True and relaunch the 
				 		# computation on the test image. If the method fails
				 		# again, "play" with the value of manual_sigma... 
				 		# and good luck.
astec_manual_sigma=15   
                        # Axial histograms fitting initialization parameter 
						# for the computation of membrane image binarization
						# axial thresholds (this parameter is used iif 
						# manual = True).
						# One may need to test different values of 
						# manual_sigma. We suggest to test values between 5 and
						# 25 in case of initial failure. Good luck.

astec_hard_thresholding=False 
                              # If the previous membrane threshold method 
							  # failed, one can force the thresholding with a
							  # "hard" threshold applied on the whole image. 
							  # To do so, this option must be set to True.
astec_hard_threshold=1.0      
                              # If hard_thresholding = True, the enhanced 
							  # membranes image is thresholded using this 
							  # parameter (value 1 seems to be ok for 
							  # time-point t001 of Aquila embryo for example).

# Tensor voting framework
astec_sigma_TV=3.6    
                      # parameter which defines the voting scale for membrane
					  # structures propagation by tensor voting method (real
					  # coordinates). 
				 	  # This parameter shoud be set between 3 um (little cells)
				 	  # and 4.5 um(big gaps in the binarized membrane image)
astec_sigma_LF=0.9    
                      # Smoothing parameter for reconstructed image (in real
					  # coordinates). It seems that the default value = 0.9 um
					  # is ok for classic use.
astec_sample=0.2      
                      # Parameter for tensor voting computation speed 
					  # optimisation (do not touch if not bewared)
astec_rayon_dil=3.6   
                      # dilatation ray for propagated ROI from time t to t+1
					  # (default: 3.6, in real coordinates) 

# Fused image conversion parameters (if astec_fusion_u8_method is set to 1)
astec_min_percentile=0.01   
                      # mc-adhocFuse parameter of type %f (default: 0.01)
astec_max_percentile=0.99   
                      # mc-adhocFuse parameter of type %f (default: 0.99)
astec_min_method='cellinterior'
                      # mc-adhocFuse param. (default: 'cellinterior')
				      # taken in global|cell|cellborder|cellinterior|voxel
astec_max_method='cellborder'
                      # mc-adhocFuse parameter (default: 'cellborder')
					  # taken in global|cell|cellborder|cellinterior|voxel
astec_sigma_hybridation=5.0 
                      # mc-adhocFuse parameter of type %f (default: 5.0)

# Default parameters (for classical use, default values should not be changed)
astec_RadiusOpening=20 		
                      # (using the approximation of a sphere of radius 20
					  # voxels as a structuring element)
astec_Thau= 25 				
                      # s(c)=h2+(c).N2(c) >t identical
astec_MinVolume=1000 		
                      # Miss Suppressing cells with to small volumes (not
					  # yet in supdata)
astec_VolumeRatioBigger=0.5 
                      # If a cell in St+1 is at least 50% bigger than its
					  # progeny in St+1, 
astec_VolumeRatioSmaller=0.1
                      # Cells in St+1 that are 10% or smaller than their
					  # equivalent in St+1 are tagged for correction
astec_MorphosnakeIterations=10 
                      # Then, an active contour algorithm is applied
				      # using the dilated shape of c, obtained by 
					  # iterating 10 times
astec_NIterations=200 		
                      # The algorithm is applied up to stability 
					  # (at th voxels) or after n iterations 
					  # (in our case th = 103 and n = 200). 
astec_DeltaVoxels=10**3  	
                      # y (at th voxels)
astec_Volum_Min_No_Seed=100 
                      # Then, if the volume of c is greater than 100 
					  # voxels (2.7 um3)
astec_nb_proc=10 			
                      # Number of processor ...
astec_nb_proc_ace=7   		
                      # number of processors for ACE (7 is recommanded)



##################################
### POST-CORRECTION PARAMETERS ###
##################################

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



