#!/usr/bin/python2.7



import os, sys, imp

assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC"))
sys.path.append(os.path.join(os.path.dirname(__file__),"ASTEC"))
assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions"))
sys.path.append(os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions"))
assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions","cpp")), "Unable to find the 'cpp' library link in %s,\
    please install properly the library and make a logical link to its bin \
    repository at the path %s."%\
    (os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions"), \
    os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions","cpp"))

from optparse import OptionParser

from nomenclature import *
from ImageHandling import imread, imsave, SpatialImage
from ASTEC.ASTEC import segmentation_propagation
from lineage import read_lineage_tree,write_lineage_tree,timeNamed,timesNamed
from lineage_test import pkl_lineage_test, imageDict


### Options parsing

parser = OptionParser()
parser.add_option("-p", "--parameters", dest="p",
                  help="python file containing parameters definition",\
                  metavar="FILE")
parser.add_option("-e", "--embryo-rep", dest="e",
                  help="path to the embryo data", metavar="PATH")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

parser.add_option("-k", "--keep-all",
                  action="store_true", dest="keepTemporaryFiles", default=False,
                  help="keep temporary file")

(options, args) = parser.parse_args()

parameters_file=options.p

print 'verbose ['+str(options.keepTemporaryFiles)+']'

### Parameters file
if not parameters_file:
    parameters_file=raw_input('Provide the parameters file: ')
try:
    assert os.path.isfile(parameters_file)
except ValueError:
    print "Provided path '%s' was not found. Exiting."%parameters_file
    sys.exit(1)

### Reading parameters from parameters_file, put in the instance p
try:
    p = imp.load_source('*', parameters_file) # p is a structure which contains
                                              # all the parameters defined in
                                              # parameters_file
except ValueError:
    print "Unable to load the source '%s'. Exiting."%parameters_file
    sys.exit(1)


### Particular options parsing: redefining p.PATH_EMBRYO

if options.e:
    if p.PATH_EMBRYO:
        test=raw_input("Caution: parameter PATH_EMBRYO '%s' defined in '%s'\
              will be overwritten by '%s'. Are you sure? (Y/n)"\
              %(p.PATH_EMBRYO,parameters_file,options.e))
        if test != 'Y' and test != 'y':
            print "Usage: option -e should be used when PATH_EMBRYO option is\
             kept empty in parameters file. Exiting."
            sys.exit(1) 
    p.PATH_EMBRYO = options.e
    # Embryo name from path: defining p.EN (also set if options.e is provided)
    if p.EN:
        test=raw_input("Caution: parameter EN '%s' defined in '%s'\
              will be overwritten. Are you sure? (Y/n)"\
              %(p.PATH_EMBRYO,parameters_file))
        if test != 'Y' and test != 'y':
            print "Usage: option '-e' should be used when PATH_EMBRYO and EN \
             options are kept empty in parameters file. Exiting."
            sys.exit(1) 
    p.EN='' # replaceFlags function will deal with empty p.EN field

### Creating some p attributes if missing (for retrocompatibility)

if not hasattr(p, 'astec_membrane_reconstruction_method'):
    p.astec_membrane_reconstruction_method = 0
if not hasattr(p, 'astec_fusion_u8_method'):
    p.astec_fusion_u8_method = 0
if not hasattr(p, 'astec_flag_hybridation'):
    p.astec_flag_hybridation = False
if not hasattr(p, 'astec_keep_reconstruct_files'):
    p.astec_keep_reconstruct_files = False
if not hasattr(p, 'astec_rayon_dil'):
    p.astec_rayon_dil = 3.6
if not hasattr(p, 'astec_nb_proc_ace'):
    p.astec_nb_proc_ace = 7

if not hasattr(p, 'astec_min_percentile'):
    p.astec_min_percentile=0.01   
if not hasattr(p, 'astec_max_percentile'):
    p.astec_max_percentile=0.99   
if not hasattr(p, 'astec_min_method'):
    p.astec_min_method='cellinterior'
if not hasattr(p, 'astec_max_method'):
    p.astec_max_method='cellborder'
if not hasattr(p, 'astec_sigma_hybridation'):
    p.astec_sigma_hybridation=5.0 

### Building paths from nomenclature.py and parameters file

path_fuse_exp = replaceFlags(path_fuse_exp, p)
print "Fused data will be searched in directory %s"%replaceFlags(path_fuse_exp,
                                                                 p)
assert os.path.isdir(path_fuse_exp), "Provided fuse directory '%s' not found"\
                                     %path_fuse_exp
path_fuse_exp_files = replaceFlags(path_fuse_exp_files, p)

path_seg = replaceFlags(path_seg, p)
path_seg_exp = replaceFlags(path_seg_exp, p)
path_seg_exp_files = replaceFlags(path_seg_exp_files, p)
path_seg_exp_lineage = replaceFlags(path_seg_exp_lineage, p)
path_seg_exp_lineage_test = replaceFlags(path_seg_exp_lineage_test, p)
path_seg_exp_reconstruct=replaceFlags(path_seg_exp_reconstruct,p)
path_seg_exp_reconstruct_files=replaceFlags(path_seg_exp_reconstruct_files,p)
if not p.astec_keep_reconstruct_files:
    path_seg_exp_reconstruct_files = None
path_log_file = replaceFlags(path_seg_logfile, p)
path_log_file = replaceEXECUTABLE(path_log_file, __file__)


dir_log_file = os.path.dirname(path_log_file)
if not os.path.isdir(dir_log_file):
        os.makedirs(dir_log_file)

### Segmentation directory and subdirectory construction

if not os.path.isdir(path_seg):
    os.mkdir(path_seg)  
if not os.path.isdir(path_seg_exp):
    os.mkdir(path_seg_exp)  
if p.astec_keep_reconstruct_files and not \
   os.path.isdir(path_seg_exp_reconstruct):
    os.mkdir(path_seg_exp_reconstruct)  


### Log file

if os.path.exists(os.getcwd()+os.path.sep+'.git'):
    os.system('echo "# astec-package version: `git describe`" >> %s'\
        %path_log_file)
else:
    os.system('echo "# astec-package version was not found, meaning it was not\
     cloned from the inria forge project" >> %s'%path_log_file)
with open(path_log_file, 'a') as log_file:
    log_file.write('# Python executable: '+sys.executable+'\n')
    log_file.write('# Working directory: %s\n'%os.getcwd())
    log_file.write('# Parameters file: %s\n'% parameters_file)
    log_file.write('# Embryo path: %s\n'% p.PATH_EMBRYO)
    log_file.write('# Embryo name: %s\n'% p.EN)
    log_file.write('# Fusion path: %s\n'% path_fuse_exp)
    log_file.write('# Command line:\n'+(' '.join(sys.argv))+'\n\n\n')

### Copy of parameters file

os.system("cp -f "+parameters_file+" "+path_seg_exp )











######################################
### Segmentation Propagation Stuff ###
######################################

# ASTEC  segmentation propagation

# Read the lineage tree (in case it was previously created)
lin_tree_information=read_lineage_tree(path_seg_exp_lineage) 

begin=p.begin+p.raw_delay
end=p.end+p.raw_delay

####SAFETY CHECK AFTER RELAUNCH
if 'lin_tree' in lin_tree_information:
    import numpy as np
    cellat={}
    for y in lin_tree_information['lin_tree']:
        t=y/10**4
        if t not in cellat:
            cellat[t]=1
        else:
            cellat[t]+=1

    restart=-1
    t=begin
    while restart==-1 and t<=end:
        time_segment=t+p.delta #Time point of Segmentation 
        segmentation_file=replaceTIME(path_seg_exp_files, time_segment) # seg
        if not os.path.isfile(segmentation_file): 
            print 'Miss segmentation file at %s -> %s'%(t, segmentation_file)
            restart=t
        else:
            if cellat[t]==0:
                    print 'Miss lineage at '+str(t)
                    restart=t
            else:
                try :
                    seg=imread(segmentation_file)
                except IOError:
                    print ' Error in '+segmentation_file
                    restart=t
        if restart==-1:
            print ' Safe segmentation at '+str(t)
        t+=1
    begin=restart
    print ' --> Restart at t%03d'%(begin)
else:
    print ' --> Start at t%03d'%(begin)


### PROCESS PROPAGATION SEGMENTATION 
for t in range(begin, end):
    time_segment=t+p.delta #Time point of Segmentation 
    print 'Starting the segmentation at ' + str(time_segment)
    fused_file_ref=os.path.join(path_fuse_exp,replaceTIME(path_fuse_exp_files, t)+".inr") #Previous image file
    fused_file=os.path.join(path_fuse_exp,replaceTIME(path_fuse_exp_files, time_segment)+".inr") #To be segmented
    segmentation_file_ref=os.path.join(path_seg_exp,replaceTIME(path_seg_exp_files, t)+".inr") #Prev. seg file
    segmentation_file=os.path.join(path_seg_exp,replaceTIME(path_seg_exp_files, time_segment)+".inr") #Output seg
    reconstruct_file=None

    if p.astec_keep_reconstruct_files:
        reconstruct_file=replaceTIME(path_seg_exp_reconstruct_files, \
                                     time_segment)
    #  TEMPORARY FOLDER
    temporary_folder=replaceTIME(os.path.join(path_seg_exp,'TEMP_'+FLAG_TIME),\
                                t)
    os.system("mkdir -p " + temporary_folder ) # Make temporary folder

    vf_file=replaceTimes( \
            os.path.join( \
            temporary_folder,'VF_t'+FLAG_TIMEREF+'_on_t'+FLAG_TIMEFLO+'.inr' \
            ), {FLAG_TIMEREF:t,FLAG_TIMEFLO:time_segment}) 
    h_min_files=replaceTIME(os.path.join(temporary_folder, \
                            'h_min_t$TIME_h$HMIN_s$SIGMA.inr'),time_segment)
    seed_file=replaceTIME(os.path.join(temporary_folder,'Seed_t$TIME.inr'),t)
    print vf_file
    print h_min_files
    print seed_file

    #PROCESS PROGATION SEGMENTATION
    seg_from_opt_h, lin_tree_information=segmentation_propagation( \
        t,fused_file_ref,segmentation_file_ref, fused_file,  \
        seed_file, vf_file, h_min_files,  \
        p.astec_h_min_min, p.astec_h_min_max, p.astec_sigma1,  \
        lin_tree_information, p.delta, p.astec_nb_proc, \
        membrane_reconstruction_method=p.astec_membrane_reconstruction_method,\
        fusion_u8_method=p.astec_fusion_u8_method, \
        flag_hybridation=p.astec_flag_hybridation, \
        RadiusOpening=p.astec_RadiusOpening, Thau=p.astec_Thau, \
        MinVolume=p.astec_MinVolume,  \
        VolumeRatioBigger=p.astec_VolumeRatioBigger, \
        VolumeRatioSmaller=p.astec_VolumeRatioSmaller, \
        MorphosnakeIterations=p.astec_MorphosnakeIterations, \
        NIterations=p.astec_NIterations, DeltaVoxels=p.astec_DeltaVoxels, \
        rayon_dil=p.astec_rayon_dil, \
        sigma_membrane=p.astec_sigma_membrane, \
        manual=p.astec_manual, \
        manual_sigma=p.astec_manual_sigma, \
        hard_thresholding=p.astec_hard_thresholding, \
        hard_threshold=p.astec_hard_threshold, \
        sensitivity=p.astec_sensitivity, \
        sigma_TV=p.astec_sigma_TV, \
        sigma_LF=p.astec_sigma_LF, \
        sample=p.astec_sample, \
        keep_membrane=False, keep_all=False,  nb_proc_ACE=p.astec_nb_proc_ace,\
        min_percentile=p.astec_min_percentile, \
        max_percentile=p.astec_max_percentile, \
        min_method=p.astec_min_method, max_method=p.astec_max_method,\
        sigma_hybridation=p.astec_sigma_hybridation, \
        path_u8_images=reconstruct_file, \
        verbose=True,
        keepTemporaryFiles=options.keepTemporaryFiles)
    
    #SAVE OUTPUT
    print 'Write the segmentation in ' + segmentation_file
    imsave(segmentation_file, seg_from_opt_h)
    #Save the current lineage tree
    write_lineage_tree(path_seg_exp_lineage,lin_tree_information)
    if options.keepTemporaryFiles==False:
        os.system("rm -rf  " + temporary_folder ) #delete temporary folder

print 'ASTEC SEGMENTATION DONE'

### PROCESS LINEAGE TREE FILE VERIFICATION
print 'PROCESS LINEAGE TREE VERIFICATION'
image_dict_seg=imageDict(path_seg_exp_files.replace(FLAG_TIME,"*"))
report=pkl_lineage_test(lin_tree_information, image_dict_seg, \
                         file_out=path_seg_exp_lineage_test)
print report
print 'LINEAGE TREE FILE VERIFICATION DONE'
