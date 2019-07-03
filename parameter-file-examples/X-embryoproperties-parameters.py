############################################################
##
## useful parameters for cell properties computation
## from a co-registered sequence
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

## ##### explanation #####
##
## time points to be processed
##
## begin: first time point
## end: last time point
##



# EXP_INTRAREG = 'RELEASE'

## ##### explanation #####
##
## Path to the sequence intra-registration result
## Fusion results are stored in <PATH_EMBRYO>/INTRAREG/INTRAREG_<EXP_INTRAREG>
## default value is
## EXP_INTRAREG = 'RELEASE'
##
## if the post-corrected images have been co-registered, they will be
## used for the cell properties computation
## (hence the image from <PATH_EMBRYO>/INTRAREG/INTRAREG_<EXP_INTRAREG>/POST/)
## else the co-registered segmentation images are used
##





######################################################################
##
## properties computation parameters (advanced)
##
######################################################################



# properties_nb_proc = None

## ##### explanation #####
##
## Number of processors used for the computation (default is 1)
## Due to the simultaneous opening of several image files, it seems
## that using parallelism for this step may be hazardous
##


