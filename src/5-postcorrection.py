#!/usr/bin/python2.7


#!/usr/bin/python2.7



import os, sys, imp

assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC"))
sys.path.append(os.path.join(os.path.dirname(__file__),"ASTEC"))
assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC", \
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
from lineage import write_tlp_from_lin_tree,read_lineage_tree,\
					write_lineage_tree,timeNamed,timesNamed
from post_correction import apply_cell_fusion,remove_too_little_branches
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

(options, args) = parser.parse_args()

parameters_file=options.p

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


### Building paths from nomenclature.py and parameters file

path_seg = replaceFlags(path_seg, p)
path_seg_exp = replaceFlags(path_seg_exp, p)

print "Segmentation data will be searched in directory %s"%replaceFlags(
															   path_seg_exp, p)
assert os.path.isdir(path_seg_exp), \
	   "Provided segmentation directory '%s' not found"%path_seg_exp

path_seg_exp_files = replaceFlags(path_seg_exp_files, p)
path_seg_exp_lineage = replaceFlags(path_seg_exp_lineage, p)

path_post = replaceFlags(path_post, p)
path_post_exp = replaceFlags(path_post_exp, p)
path_post_exp_files = replaceFlags(path_post_exp_files, p)
path_post_exp_lineage = replaceFlags(path_post_exp_lineage, p)
path_post_exp_lineage_test = replaceFlags(path_post_exp_lineage_test, p)
path_log_file = replaceFlags(path_post_logfile, p)

### Segmentation directory and subdirectory construction

if not os.path.isdir(path_post):
    os.mkdir(path_post)  
if not os.path.isdir(path_post_exp):
    os.mkdir(path_post_exp)  


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
    log_file.write('# Segmentation path: %s\n'% path_seg_exp)
    log_file.write('# Command line:\n'+(' '.join(sys.argv))+'\n\n\n')

### Copy of parameters file

os.system("cp -f "+parameters_file+" "+path_post_exp )







##########################################
### Segmentation Post-Correction Stuff ###
##########################################

# Post correction of the segmentation

# Read the lineage tree (in case it was previously created)
lin_tree_information=read_lineage_tree(path_seg_exp_lineage) 



# Lineage Tree Post-correction

### CORRECT THE LINEAGE TREE WITH THE CELLS VOLUMES
lin_tree_cor, new_volumes, to_fuse, been_fused=\
	remove_too_little_branches(lin_tree_information['lin_tree'], \
		lin_tree_information['volumes_information'], \
		p.postcor_Volume_Threshold, soon=p.postcor_Soon)


### APPLYING THE CORRECTION ON THE IMAGES
apply_cell_fusion(lin_tree_information['lin_tree'], \
	   lin_tree_information['volumes_information'], \
	   to_fuse,os.path.join(path_seg_exp,path_seg_exp_files), \
                  os.path.join(path_post_exp,path_post_exp_files),p.begin+p.raw_delay, \
	   p.end+p.raw_delay, p.delta)

#SAVE THE NEW LINEAGE TREE
lin_tree_information['lin_tree']=lin_tree_cor


### SAVE THE FINAL LINEAGE TREE
write_lineage_tree(path_post_exp_lineage,lin_tree_information)

### PROCESS LINEAGE TREE FILE VERIFICATION
print 'PROCESS LINEAGE TREE VERIFICATION'
image_dict_seg_post=imageDict(path_post_exp_files.replace("$TIME","*"))
report=pkl_lineage_test(lin_tree_information, image_dict_seg_post, \
						file_out=path_post_exp_lineage_test)
print report
print 'LINEAGE TREE VERIFICATION DONE'

