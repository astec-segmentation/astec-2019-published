############################################################
##
## useful parameters for sequence intra-registration
##
## to be used with 1.5-intraregistration.py
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

## ##### explanation #####
##
## time points to be processed
##
## begin: first time point
## end: last time point
## delta: delta/time interval between two time points (if one does not want
##        to deal with every single time point) (default = 1)
##



# EXP_FUSE = 'RELEASE'

## ##### explanation #####
##
## Path to the fusion result
## Fusion results are stored in <PATH_EMBRYO>/FUSE/FUSE_<EXP_FUSE>
## default value is
## EXP_FUSE = 'RELEASE'
##



# EXP_SEG = 'RELEASE'

## ##### explanation #####
##
## Path to the segmentation result
## n results are stored in <PATH_EMBRYO>/SEG/SEG_<EXP_SEG>
## default value is
## EXP_SEG = 'RELEASE'
##
## used when "intra_registration_template_type = 'SEGMENTATION'"
## 



# EXP_POST = 'RELEASE'

## ##### explanation #####
##
## Path to the post-segmentation result
## n results are stored in <PATH_EMBRYO>/POST/POST_<EXP_POST>
## default value is
## EXP_POST = 'RELEASE'
##
## used when "intra_registration_template_type = 'POST-SEGMENTATION'"
## 



# EXP_INTRAREG = 'RELEASE'

## ##### explanation #####
##
## Path to the sequence intra-registration result
## Fusion results are stored in <PATH_EMBRYO>/INTRAREG/INTRAREG_<EXP_INTRAREG>
## default value is
## EXP_INTRAREG = 'RELEASE'



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
## intra-registration parameters (advanced)
##
######################################################################



# intra_registration_compute_registration = True
# intra_registration_transformation_type = 'rigid'
# intra_registration_transformation_estimation_type = 'wlts'
# intra_registration_lts_fraction = 0.55
# intra_registration_pyramid_highest_level = 6
# intra_registration_pyramid_lowest_level = 3
# intra_registration_normalization = True

## ##### explanation #####
##
## parameters for the co-registration of successive fused images
##



# intra_registration_resolution = 0.6
# intra_registration_reference_index = None
# intra_registration_template_type = "FUSION"
# intra_registration_template_threshold = None
# intra_registration_margin = None

## ##### explanation #####
##
## parameters for template building
##
## intra_registration_resolution: gives the resulting (isotropic) voxel size
##   (as the 'target_resolution' gives the voxel size of the fused images). However, for
##   visualization purposes, it may be indicated to have a larger voxel size (hence the 0.6
##   instead the 0.3)
## intra_registration_reference_index: defines the still image after drift compensation
## intra_registration_template_type:  'FUSION' | 'SEGMENTATION' | 'POST-SEGMENTATION'
##   The template is built so that the useful information of all resampled images fits into it.
##   Useful information can be issued from either the fused sequence, the segmentation sequence or
##   the post-segmentation sequence. 
## intra_registration_template_threshold:
##   Giving a threshold with the 'intra_registration_template_type', only points above the 
##   threshold are considered to be included in the template after resampling, this allows 
##   to reduce the template. According the background value is either 0 or 1 in both the 
##   segmentation and the post-segmentation sequences, setting this threshold to 2 for these 
##   sequences allows to keep the entire embryo in the resampled/reconstructed sequence.
## intra_registration_margin:
##   In addition, a margin can be given for a more comfortable visualization. By default, it is
##   0 when only fusion images are used, and 10 if either segmentation or post-segmentation 
##   images are also used
##



# intra_registration_resample_fusion_images = True
# intra_registration_resample_segmentation_images = False
# intra_registration_resample_post_segmentation_images = False

# intra_registration_sigma_segmentation_images = 1.0

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

## ##### explanation #####
##
## outputs
## 
## intra_registration_resample_fusion_images:
## intra_registration_resample_segmentation_images:
## intra_registration_resample_post_segmentation_images:
##   indicates whether the [fusion,segmentation,post-segmentation] images
##   have to be resampled in the template geometry
##   If required, it assumes that segmentation and post-segmentation images exists
##   and are respectively located in <PATH_EMBRYO>/SEG/SEG_<EXP_SEG> and
##   <PATH_EMBRYO>/POST/POST_<EXP_POST>
## intra_registration_sigma_segmentation_images:
##   gives the sigma (gaussian standard deviation) used to resample the
##   cell/label images (segmentation and/or post-segmentation)
## intra_registration_[xy,xz,yz]_movie_fusion_images:
## intra_registration_[xy,xz,yz]_movie_segmentation_images:
## intra_registration_[xy,xz,yz]_movie_post_segmentation_images:
##   indicates the index(es) of the [xy,xz,yz] sections to build
##   [xy,xz,yz]+t movies
##

