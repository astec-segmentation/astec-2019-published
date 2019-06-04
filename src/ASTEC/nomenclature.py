#
# File defining the nomenclature for ASTEC experiments
#

import os


#
#
# General definition
#
#

FLAG_EXECUTABLE = '$EXECUTABLE'

#
# main path = path to embryo directory
# - rawdata have to be a sub-directory of it
# - results will be sub-directories also
#

FLAG_PATH_EMBRYO = '$PATHEMBRYO'

#
# Embryo Name
# CRBM convention is: YYMMDD-SaintOfTheDays-Stage
#
FLAG_EN = '$EN'


#
# Definition of main sub-directories of the embryo directory
# apart the directory containing the raw data
# they are created by the ASTEC stages/utilities
#
# - RAW DATA are stored in [FLAG_PATH_EMBRYO]/[FLAG_DIR_RAWDATA]
# - Fusion experiments stored in
#   [FLAG_PATH_EMBRYO]/[DIR_STAGE_FUSE]/[DIR_STAGE_FUSE]_[FLAG_EXP_FUSE]
# - Intra-registration experiments stored in
#   [FLAG_PATH_EMBRYO]/[DIR_STAGE_INTRAREG]/[DIR_STAGE_INTRAREG]_[FLAG_EXP_FUSE]
# - Mars experiments (segmentation of first time point) stored in
#   [FLAG_PATH_EMBRYO]/[DIR_STAGE_SEG]/[DIR_STAGE_SEG]_[FLAG_EXP_SEG]
# - Segmentation propagation experiments stored in
#   [FLAG_PATH_EMBRYO]/[DIR_STAGE_SEG]/[DIR_STAGE_SEG]_[FLAG_EXP_SEG]
# - Post correction experiments stored in
#   [FLAG_PATH_EMBRYO]/[DIR_STAGE_SEG]/[DIR_STAGE_SEG]_[FLAG_EXP_POST]

FLAG_DIR_RAWDATA = '$DIRRAWDATA'
DIR_STAGE_FUSE = 'FUSE'
DIR_STAGE_INTRAREG = 'INTRAREG'
DIR_STAGE_SEG = 'SEG'
DIR_STAGE_POST = 'POST'

DIR_STAGE_REG = 'REG'

FLAG_EXP_FUSE = '$FUSE'		# defining fusion experiments subdirectory
FLAG_EXP_INTRAREG = '$INTRAREG'			# defining registration experiments subdirectory
FLAG_EXP_SEG = '$SEG'			# defining seg propagation experiments subdirectory
FLAG_EXP_POST = '$POST'		# defining post correction experiments subdirectory

FLAG_EXP_REG = '$REG'

#
# Flags relative to file name
#

FLAG_TIMESTAMP = '$TIMESTAMP'  # time stamp of the experiment


FLAG_TIME = '$TIME'			# Time-point of an embryo snapshot
FLAG_TIMEREF = '$TIMEREF'		# For registration,time-point of reference snapshot
FLAG_TIMEFLO = '$TIMEFLO'		# For registration,time-point of floating snapshot

FLAG_RESULT_IMAGE_SUFFIX_FUSE = '$RESULTIMAGESUFFIX'


#
# raw data specific definitions: FLAGS
#
# definitions of sub-directories in [PATHEMBRYO]/[FLAG_DIR_RAWDATA]
#
# if not defined, they will be
# LC/Stack0000
# LC/Stack0001
# RC/Stack0000
# RC/Stack0001
#
# else
# - LC/Stack0000 images are in PATHEMBRYO]/[FLAG_DIR_RAWDATA]/[FLAG_DIR_LEFTCAM_STACKZERO]
# - LC/Stack0001 images are in PATHEMBRYO]/[FLAG_DIR_RAWDATA]/[FLAG_DIR_LEFTCAM_STACKONE]
# - RC/Stack0000 images are in PATHEMBRYO]/[FLAG_DIR_RAWDATA]/[FLAG_DIR_RIGHTCAM_STACKZERO]
# - RC/Stack0001 images are in PATHEMBRYO]/[FLAG_DIR_RAWDATA]/[FLAG_DIR_RIGHTCAM_STACKone]
#
# Images names are defined by
# [FLAG_NAME_LEFT_RAW_IMAGE]${TIME} and [FLAG_NAME_RIGHT_RAW_IMAGE]${TIME}
# default are 'Cam_Left_00' and 'Cam_Right_00' (time is assumed to be encoded on 3 bytes)
#

FLAG_DIR_LEFTCAM_STACKZERO = '$DIRLEFTCAMSTACKZERO'
FLAG_DIR_LEFTCAM_STACKONE = '$DIRLEFTCAMSTACKONE'
FLAG_DIR_RIGHTCAM_STACKZERO = '$DIRFIGHTCAMSTACKZERO'
FLAG_DIR_RIGHTCAM_STACKONE = '$DIRFIGHTCAMSTACKONE'

FLAG_NAME_LEFT_RAW_IMAGE = '$NAMELEFTRAWIMAGE'
FLAG_NAME_RIGHT_RAW_IMAGE = '$NAMERIGHTRAWIMAGE'

#
# raw data specific definitions: variables
#
# angle #1 = left camera, stack 0
# angle #2 = right camera, stack 0
# angle #3 = left camera, stack 1
# angle #4 = right camera, stack 1
#

path_rawdata = os.path.join(FLAG_PATH_EMBRYO, FLAG_DIR_RAWDATA)

# 1st image from the left camera
path_rawdata_angle1 = os.path.join(path_rawdata, FLAG_DIR_LEFTCAM_STACKZERO)

# 1st image from the right camera
path_rawdata_angle2 = os.path.join(path_rawdata, FLAG_DIR_RIGHTCAM_STACKZERO)

# 2nd image from the left camera
path_rawdata_angle3 = os.path.join(path_rawdata, FLAG_DIR_LEFTCAM_STACKONE)

# 2nd from the right camera
path_rawdata_angle4 = os.path.join(path_rawdata, FLAG_DIR_RIGHTCAM_STACKONE)


path_rawdata_angle1_files = FLAG_NAME_LEFT_RAW_IMAGE + FLAG_TIME
path_rawdata_angle2_files = FLAG_NAME_RIGHT_RAW_IMAGE + FLAG_TIME
path_rawdata_angle3_files = FLAG_NAME_LEFT_RAW_IMAGE + FLAG_TIME
path_rawdata_angle4_files = FLAG_NAME_RIGHT_RAW_IMAGE + FLAG_TIME


#
# FUSION DATA
#


# path fused images
path_fuse = os.path.join(FLAG_PATH_EMBRYO, DIR_STAGE_FUSE)
# path for fusion workspace
path_fuse_exp = os.path.join(path_fuse, DIR_STAGE_FUSE + '_' + FLAG_EXP_FUSE)
# fused images generic name
path_fuse_exp_files = FLAG_EN + '_fuse_t' + FLAG_TIME

# log files
path_fuse_logdir = os.path.join(path_fuse_exp, 'LOGS')
path_fuse_historyfile = os.path.join(path_fuse_logdir, FLAG_EXECUTABLE + '-history.log')
path_fuse_logfile = os.path.join(path_fuse_logdir, FLAG_EXECUTABLE + '-' + FLAG_TIMESTAMP + '.log')

#
# INTRA REGISTRATION DATA
#

# main path
path_intrareg = os.path.join(FLAG_PATH_EMBRYO, DIR_STAGE_INTRAREG)
# path for intra-registration workspace
path_intrareg_exp = os.path.join(path_intrareg, DIR_STAGE_INTRAREG + '_' + FLAG_EXP_INTRAREG)

# path for intra-registration co-registration transformation
path_intrareg_cotrsf = os.path.join(path_intrareg_exp, 'CO-TRSFS')
path_intrareg_cotrsf_files = FLAG_EN + '_intrareg_flo' + FLAG_TIMEFLO + '_ref' + FLAG_TIMEREF + '.trsf'
path_intrareg_trsf = os.path.join(path_intrareg_exp, 'TRSFS')
path_intrareg_trsf_files = FLAG_EN + '_intrareg_t' + FLAG_TIME + '.trsf'

path_intrareg_fuse = os.path.join(path_intrareg_exp, 'FUSE')
path_intrareg_seg = os.path.join(path_intrareg_exp, 'SEG')
path_intrareg_post = os.path.join(path_intrareg_exp, 'POST')

# log files
path_intrareg_logdir = os.path.join(path_intrareg_exp, 'LOGS')
path_intrareg_historyfile = os.path.join(path_intrareg_logdir, FLAG_EXECUTABLE + '-history.log')
path_intrareg_logfile = os.path.join(path_intrareg_logdir, FLAG_EXECUTABLE + '-' + FLAG_TIMESTAMP + '.log')

# resampled images
path_intrareg_fuse = os.path.join(path_intrareg_exp, 'FUSE')
path_intrareg_fuse_files = FLAG_EN + '_intrareg_fuse_t' + FLAG_TIME

path_intrareg_seg = os.path.join(path_intrareg_exp, 'SEG')
path_intrareg_seg_files = FLAG_EN + '_intrareg_seg_t' + FLAG_TIME

path_intrareg_post = os.path.join(path_intrareg_exp, 'POST')
path_intrareg_post_files = FLAG_EN + '_intrareg_post_t' + FLAG_TIME

# movies
path_intrareg_movies = os.path.join(path_intrareg_exp, 'MOVIES')


#
# MARS DATA
#


# path mars images
# path_mars = os.path.join(FLAG_PATH_EMBRYO, DIR_STAGE_MARS)
# path for mars workspace
# path_mars_exp = os.path.join(path_mars, DIR_STAGE_MARS + '_' + FLAG_EXP_MARS)
# mars images names
path_mars_exp_files = FLAG_EN + '_mars_t' + FLAG_TIME

# path_mars_exp_reconstruct = os.path.join(path_mars_exp, "RECONSTRUCTION")
# mars temporary images file names (for reconstruction)
# path_mars_exp_reconstruct_files = os.path.join(path_mars_exp_reconstruct, FLAG_EN + '_rec_t$TIME.inr')
# logfile
# path_mars_logdir = os.path.join(path_mars_exp, 'LOGS')
# path_mars_historyfile = os.path.join(path_mars_logdir, FLAG_EXECUTABLE + '-history.log')
# path_mars_logfile = os.path.join(path_mars_logdir, FLAG_EXECUTABLE + '-' + FLAG_TIMESTAMP + '.log')


#
# SEGMENTATION DATA
#


# path segmented images
path_seg = os.path.join(FLAG_PATH_EMBRYO, DIR_STAGE_SEG) 
# path for segmentation workspace
path_seg_exp = os.path.join(path_seg, DIR_STAGE_SEG + '_' + FLAG_EXP_SEG)
# segmentated images names
# path_seg_exp_files = os.path.join(path_seg_exp, FLAG_EN + '_seg_t$TIME.inr')
path_seg_exp_files = FLAG_EN + '_seg_t' + FLAG_TIME

# log files
path_seg_logdir = os.path.join(path_seg_exp, 'LOGS')
path_seg_historyfile = os.path.join(path_seg_logdir, FLAG_EXECUTABLE + '-history.log')
path_seg_logfile = os.path.join(path_seg_logdir, FLAG_EXECUTABLE + '-' + FLAG_TIMESTAMP + '.log')


# cells lineage file
# path_seg_exp_lineage = os.path.join(path_seg_exp, FLAG_EN + '_seg_lineage.pkl')
path_seg_exp_lineage = FLAG_EN + '_seg_lineage.pkl'

# segmentation temporary images path (for reconstruction)
path_seg_exp_reconstruct = os.path.join(path_seg_exp, "RECONSTRUCTION")
# segmentation temporary images file name (for reconstruction)
path_seg_exp_reconstruct_files = os.path.join(path_seg_exp_reconstruct, FLAG_EN + '_rec_t$TIME.inr')
# cells lineage test file
path_seg_exp_lineage_test = os.path.join(path_seg_exp, FLAG_EN + '_seg_lineage.test')
# logfiles
# path_seg_logfile = os.path.join(path_seg_exp, '4-astec.log')
# logfile


#
# MANUAL CORRECTION DATA
#


# logfiles
path_mancor_logfile = os.path.join(path_seg_exp, '3-manualcorrection.log')

#
# POST CORRECTION DATA
#


# path post corrected images
path_post = os.path.join(FLAG_PATH_EMBRYO, DIR_STAGE_POST)
# path for post correction workspace
path_post_exp = os.path.join(path_post, DIR_STAGE_POST + '_' + FLAG_EXP_POST)
# post corrected images names
path_post_exp_files = FLAG_EN + '_post_t' + FLAG_TIME
# path_post_exp_files = os.path.join(path_post_exp, FLAG_EN + '_post_t$TIME.inr')

# post corrected cells lineage file
path_post_exp_lineage = os.path.join(path_post_exp, FLAG_EN + '_post_lineage.pkl')
# post corrected cells lineage test file
path_post_exp_lineage_test = os.path.join(path_post_exp, FLAG_EN + '_post_lineage.test')
# logfile
path_post_logfile = os.path.join(path_post_exp, '5-postcorrection.log')


#
# INTRA REGISTRATION DATA
#


#path_intrareg= os.path.join(FLAG_PATH_EMBRYO, DIR_STAGE_REG) # Path intra registration data
#path_intrareg_exp= os.path.join(path_intrareg, DIR_STAGE_REG+'_'+FLAG_EXP_REG) # Path intra registration data
#path_intrareg_step_files= os.path.join(path_intrareg_exp, \
#		FLAG_EN+'_reg_t$TIMEFLO_t$TIMEREF.trsf') # Intra registration step-by-step
											# trsf file names
#path_intrareg_multiple_files= os.path.join(path_intrareg, \
#		FLAG_EN+'_reg_compose_t$TIME_t$TIME.trsf') # Intra registration composed
											  #trsf file names
#path_intrareg_change_files= os.path.join(path_intrareg, \
#		FLAG_EN+'_reg_compose_t$TIME.trsf') # Intra registration recentered trsf
									   # file names
#path_intrareg_change_template= os.path.join(path_intrareg, \
#		FLAG_EN+'_reg_compose_template.inr.gz') # Intra registration template file
										   # name for recentered trsfs
#iso_intra_registration=1.0		# Parameter for intra registration resampled
								# images resolution
















#postsegment_files=datapath+"GLACE/SEG/POST/"+EN+'_glas_seg_post_t$TIME.inr' #Segmentation output files

#INTRA REGISTRATION COMPOSED WITH ROTATION SO THAT GERMINAL CELLS ARE AT THE DOWN OF THE IMAGE # NOT USED NOW
intrareg_germinal_Path=path_intrareg_exp+"COMPOSE_GERMINAL/" # path intra registration data under germinal cells reorientation constraint
intrareg_germinal_file=intrareg_germinal_Path+FLAG_EN+'_germinal.trsf' # transformation (rotation) to compose with all the intra-registration trsf files 
intrareg_germinal_files=intrareg_germinal_Path+FLAG_EN+'_germinal_t$TIME.trsf' #  intra registration recentered trsf file names
intrareg_germinal_template=intrareg_germinal_Path+FLAG_EN+'_germinal_template.inr.gz' #  intra registration template file name for recentered trsfs




#RECONSTRUCTION DATA 
reconstruct_Path= os.path.join(path_fuse_exp,"RECONSTRUCTION") #path reconstructed images
reconstruct_files= os.path.join(reconstruct_Path, FLAG_EN+'_rec_t$TIME.inr') #  reconstructed images names



#SEGMENTED DATA 
segmented_Path= os.path.join(FLAG_PATH_EMBRYO, DIR_STAGE_SEG) #segmented images
mars_file=segmented_Path+FLAG_EN+'_fuse_mars_t$TIME.inr' #Segmentation output files
segmentation_files=segmented_Path+FLAG_EN+'_fuse_seg_t$TIME.inr' #Segmentation output files
lineage_tree_filename=segmented_Path+FLAG_EN+'_fuse_seg_lin_tree.pkl' #The main lineage tree file output 
lineage_tree_test_filename=segmented_Path+FLAG_EN+'_fuse_seg_lin_tree.test' #The main lineage tree test file output 

#MAPPING FILE
mapping_path=segmented_Path+FLAG_EN+"_fuse_seg_t$TIME.map"

#POST SEGMENTED DATA 
postsegment_Path=segmented_Path+"POST/" #post segmentation images
postsegment_files=postsegment_Path+FLAG_EN+'_fuse_seg_post_t$TIME.inr' #Segmentation output files
post_lineage_tree_filename=postsegment_Path+FLAG_EN+'_fuse_seg_post_lin_tree.pkl' #The main lineage tree file output 
post_lineage_tree_test_filename=postsegment_Path+FLAG_EN+'_fuse_seg_post_lin_tree.test' #The main lineage tree test file output


#CELL NAMES
name_file= os.path.join(FLAG_PATH_EMBRYO, FLAG_EN+"-names.txt")
named_lineage_tree_filename=postsegment_Path+FLAG_EN+'_fuse_seg_post_lin_tree_named.pkl' #The main lineage tree file output 
named_lineage_tree_test_filename=postsegment_Path+FLAG_EN+'_fuse_seg_post_lin_tree_named.test' #The main lineage tree test file output 


###################################################################
####### VIRTUAL EMBRYO SHOULD BE TAKEN OUT OF ASTEC PACKAGE #######
###################################################################

#VIRTUAL EMBRYO
#vo_login="emmanuel.faure"
#vo_passwd="ascidie"
#vo_lib_Path=astec_Path+"3DCloudEmbryo/" 		   # Path for scripts library
#vo_Path=postsegment_Path+"VO/" 					# Output path for mesh data
#vo_files=vo_Path+EN+'_fuse_seg_post_vo_t$TIME.obj' # Output files





##############################################################
################ FUNCTIONS FOR FLAGS REPLACEMENT #############
##############################################################


def replaceTIME(filename, time, flag=FLAG_TIME, timeDigits=3):
    """
    Replaces all the occurrences of "$TIME" by its value given by time of 
    type int at the format specified by _format argument (default="%03d")
    """
    time_point= '{:0{width}d}'.format(time, width=timeDigits)
    return filename.replace(flag, time_point)


#def replaceTIMEFLO(filename,time,_format="%03d"):
#    """
#    Replaces all the occurrences of "$TIMEFLO" by its value given by time of
#    type int at the format specified by _format argument (default="%03d")
#    """
#    time_point=_format%time #Format time on 3 digit
#    return filename.replace(FLAG_TIMEFLO, time_point)

#def replaceTIMEREF(filename,time,_format="%03d"):
#    """
#    Replaces all the occurrences of "$TIMEREF" by its value given by time of
#    type int at the format specified by _format argument (default="%03d")
#    """
#    time_point=_format%time #Format time on 3 digit
#    return filename.replace(FLAG_TIMEREF, time_point)


def replaceTimes(filename, d, _format="%03d"):
    """
    Replaces all the occurrences of each key from the dictionnary d given in 
    parameter with its corresponding value of type int at the format 
    specified by _format argument (default="%03d")
    """
    for k,time in d.iteritems():
        assert type(k)==str, "Unexpected non str dictionnary key '%s'"%str(k)
        if type(time)!=int:
            print "Non int val '%s' specified for key '%s'. Trying a cast."\
            	  %(str(time),k)
            time=int(time)
        assert filename.find(k)>=0
        time_point=_format%time
        filename=filename.replace(k, time_point)
    return filename


def replaceEN(filename,embryoname, check_name=False):
    """
    Replaces all the occurrences of "$EN" by its value given by EN of type str
    check_name: parameter enabling to check the EN name composition wrt the
    			CRBM nomenclature (default=False)
    """
    if check_name:
        assert embryoname.count('-')!=2, 'ERROR in embryo name %s\n\
            -> EN should follow the template YYMMDD-SaintOfTheDays-Stage'\
            %embryoname
    return filename.replace(FLAG_EN, embryoname)


def replacePATH(filename, path):
    """
    Replaces all the occurrences of "$EMBRYOPATH" by its value given by path of
    type str
    """
    #assert path, "Specified embryo path should not be empty."
    return filename.replace(FLAG_PATH_EMBRYO, path.rstrip(os.path.sep))


def replaceTIMESTAMP(filename, timestamp):
    """
    Replaces all the occurrences of "$TIMESTAMP" by its value given by EN of
    type str
    """
    import time
    if timestamp == None:
        timestamp=time.localtime()
    d=time.strftime("%Y-%m-%d-%H:%M:%S", timestamp)
    return filename.replace(FLAG_TIMESTAMP, d)


def replaceEXECUTABLE(filename, executable):
    """
    Replaces all the occurrences of "$EXECUTABLE" by its value given by EN of
    type str
    """
    import time
    if executable == None:
        local_executable = 'unknown'
    else:
        local_executable = os.path.basename(executable)
        if local_executable[-3:] == '.py':
            local_executable = local_executable[:-3]
    return filename.replace(FLAG_EXECUTABLE, local_executable)


def _genericReplaceFlag( name, flag, parameters, attributeName, defaultValue):
    attributeValue=getattr( parameters, attributeName, defaultValue)
    return name.replace(flag, attributeValue)


#
#
#


def replaceFlags(filename, parameters, timestamp=None, check_name=False):
    """
    Function that replaces in the specified str filename the flags found
    (which verify the regular expression r'\$[A-Z]*', e.g. '$EN' or $EXP_SEG)
    with the corresponding str found in the parameters (of type 'module'):
    	FLAG_PATH_EMBRYO -> parameters.PATH_EMBRYO
    	FLAG_EN -> parameters.EN
    	FLAG_EXP_FUSE -> parameters.EXP_FUSE
    	FLAG_EXP_REG -> parameters.EXP_REG
    	FLAG_EXP_SEG -> parameters.EXP_SEG
    	FLAG_EXP_POST -> parameters.EXP_POST
    """
    proc=replaceFlags
    import re
    found=re.findall(r'\$[A-Z]*',filename)

    for flag in found:

        if flag == FLAG_PATH_EMBRYO:
            if hasattr(parameters, 'PATH_EMBRYO'):
                embryo_path = parameters.PATH_EMBRYO
            else:
                embryo_path = os.getcwd()
            filename=replacePATH(filename, embryo_path)

        elif flag == FLAG_EN:
            if hasattr(parameters, 'EN'):
                embryo_name = parameters.EN
            else:
                if hasattr(parameters, 'PATH_EMBRYO'):
                    embryo_path = parameters.PATH_EMBRYO
                else:
                    embryo_path = os.getcwd()
                embryo_name = embryo_path.split(os.path.sep)[-1]
            filename = replaceEN( filename, embryo_name, check_name=check_name)

        elif flag == FLAG_DIR_RAWDATA:
            filename = _genericReplaceFlag(filename, flag, parameters, 'DIR_RAWDATA', "RAWDATA")

        elif flag == FLAG_DIR_LEFTCAM_STACKZERO:
            d = "LC" + os.path.sep + "Stack0000"
            filename = _genericReplaceFlag(filename, flag, parameters, 'DIR_LEFTCAM_STACKZERO', d)

        elif flag == FLAG_DIR_LEFTCAM_STACKONE:
            d = "LC" + os.path.sep + "Stack0001"
            filename = _genericReplaceFlag(filename, flag, parameters, 'DIR_LEFTCAM_STACKONE', d)

        elif flag == FLAG_DIR_RIGHTCAM_STACKZERO:
            d = "RC" + os.path.sep + "Stack0000"
            filename = _genericReplaceFlag(filename, flag, parameters, 'DIR_RIGHTCAM_STACKZERO', d)

        elif flag == FLAG_DIR_RIGHTCAM_STACKONE:
            d = "RC" + os.path.sep + "Stack0001"
            filename = _genericReplaceFlag(filename, flag, parameters, 'DIR_RIGHTCAM_STACKONE', d)


        elif flag == FLAG_NAME_LEFT_RAW_IMAGE:
            filename = _genericReplaceFlag(filename, flag, parameters, 'acquisition_leftcam_image_prefix', 'Cam_Left_00')

        elif flag == FLAG_NAME_RIGHT_RAW_IMAGE:
            filename = _genericReplaceFlag(filename, flag, parameters, 'acquisition_rightcam_image_prefix', 'Cam_Right_00')


        elif flag == FLAG_EXP_FUSE:
            filename = _genericReplaceFlag(filename, flag, parameters, 'EXP_FUSE', 'RELEASE')

        elif flag == FLAG_EXP_REG:
            filename = _genericReplaceFlag(filename, flag, parameters, 'EXP_REG', 'RELEASE')

        elif flag == FLAG_EXP_INTRAREG:
            filename = _genericReplaceFlag(filename, flag, parameters, 'EXP_INTRAREG', 'RELEASE')

        elif flag == FLAG_EXP_SEG:
            filename = _genericReplaceFlag(filename, flag, parameters, 'EXP_SEG', 'RELEASE')

        elif flag == FLAG_EXP_POST:
            filename = _genericReplaceFlag(filename, flag, parameters, 'EXP_POST', 'RELEASE')

        elif flag == FLAG_RESULT_IMAGE_SUFFIX_FUSE:
            filename = _genericReplaceFlag(filename, flag, parameters, 'RESULT_IMAGE_SUFFIX_FUSE', 'inr')

        elif flag == FLAG_TIMESTAMP:
            filename = replaceTIMESTAMP(filename, timestamp)

        elif flag == FLAG_TIME or flag == FLAG_TIMEFLO or flag == FLAG_TIMEREF or flag == FLAG_EXECUTABLE:
            #
            # just not to get the default message
            # $TIME may be replaced later
            #
            pass


        else:
            print "replaceFlags: flag '" + str(flag) + "' was not replaced"

    return filename