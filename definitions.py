# File defining the nomenclature

#FIND PATH AND EMBRYO NAME
import sys,os

astec_Path=os.getcwd()
#astec_Path="/home/gmicheli/TEST_ASTEC/ASTEC/160708-Aquila-St8/ASTEC-180326"
astec_Path="/home/gmicheli/MyData/DIGEM/TEST_ASTEC/RAW/171107-Karine-St8/astec-package"

tab_Path=astec_Path.split('/')
EN=tab_Path[len(tab_Path)-2]  #Embryo Name
if  EN.count('-')!=2:
	print 'ERROR in astec path '+astec_Path+ ' or embryo name '+EN
	print  '-> path should be /media/DATA/<embryoname>/ASTEC-YYMMDD'
	print '    where embryoname=YYMMDD-SaintOfTheDays-Stage'
	quit()

#datapath=astec_Path[:astec_Path.rfind('/')+1]
#astec_Path+='/'

print 'Embryo '+EN+' from '+str(begin)+' to '+str(end)

# ASTEC LIBARRIES
sys.path.append(astec_Path+"ASTEC") #Add the ASTEC Function
sys.path.append(astec_Path+'ASTEC/CommunFunctions')



#Image Path definition
rawdata_Path=datapath+"RAWDATA/"

path_angle1=rawdata_Path+"LC/Stack0000" # 1st image from the left camera, good quality image at the beginning
path_angle2=rawdata_Path+"RC/Stack0000" # 1st image from the right camera
path_angle3=rawdata_Path+"LC/Stack0001" # 2nd image from the left camera
path_angle4=rawdata_Path+"RC/Stack0001" # 2nd from the right camera


#FUSION DATA 
fuse_Path=datapath+"FUSE/" #path fused images
fused_files=fuse_Path+EN+'_fuse_t$TIME.inr' #  fused images names

#INTRA REGISTRATION DATA 
intrareg_Path=datapath+"FUSE/REG/" #path intra registration data
intrareg_step_files=intrareg_Path+EN+'_reg_t$TIMEFLO_t$TIMEREF.trsf' #  intra registration step-by-step trsf file names
intrareg_multiple_files=intrareg_Path+EN+'_reg_compose_t$TIME_t$TIME.trsf' #  intra registration composed trsf file names
intrareg_change_files=intrareg_Path+EN+'_reg_compose_t$TIME.trsf' #  intra registration recentered trsf file names
intrareg_change_template=intrareg_Path+EN+'_reg_compose_template.inr.gz' #  intra registration template file name for recentered trsfs
iso_intra_registration=1.0

#postsegment_files=datapath+"GLACE/SEG/POST/"+EN+'_glas_seg_post_t$TIME.inr' #Segmentation output files

#INTRA REGISTRATION COMPOSED WITH ROTATION SO THAT GERMINAL CELLS ARE AT THE DOWN OF THE IMAGE # NOT USED NOW
intrareg_germinal_Path=intrareg_Path+"COMPOSE_GERMINAL/" # path intra registration data under germinal cells reorientation constraint
intrareg_germinal_file=intrareg_germinal_Path+EN+'_germinal.trsf' # transformation (rotation) to compose with all the intra-registration trsf files 
intrareg_germinal_files=intrareg_germinal_Path+EN+'_germinal_t$TIME.trsf' #  intra registration recentered trsf file names
intrareg_germinal_template=intrareg_germinal_Path+EN+'_germinal_template.inr.gz' #  intra registration template file name for recentered trsfs

#RECONSTRUCTION DATA 
reconstruct_Path=fuse_Path+"RECONSTRUCTION/" #path reconstructed images
reconstruct_files=reconstruct_Path+EN+'_rec_t$TIME.inr' #  reconstructed images names

#SEGMENTED DATA 
segmented_Path=fuse_Path+"SEG/" #segmented images
mars_file=segmented_Path+EN+'_fuse_mars_t$TIME.inr' #Segmentation output files
segmentation_files=segmented_Path+EN+'_fuse_seg_t$TIME.inr' #Segmentation output files
lineage_tree_filename=segmented_Path+EN+'_fuse_seg_lin_tree.pkl' #The main lineage tree file output 
lineage_tree_test_filename=segmented_Path+EN+'_fuse_seg_lin_tree.test' #The main lineage tree test file output 

#MAPPING FILE
mapping_path=segmented_Path+EN+"_fuse_seg_t$TIME.map"

#POST SEGMENTED DATA 
postsegment_Path=segmented_Path+"POST/" #post segmentation images
postsegment_files=postsegment_Path+EN+'_fuse_seg_post_t$TIME.inr' #Segmentation output files
post_lineage_tree_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree.pkl' #The main lineage tree file output 
post_lineage_tree_test_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree.test' #The main lineage tree test file output


#CELL NAMES
name_file=datapath+EN+"-names.txt"
named_lineage_tree_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree_named.pkl' #The main lineage tree file output 
named_lineage_tree_test_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree_named.test' #The main lineage tree test file output 


###################################################################
####### VIRTUAL EMBRYO SHOULD BE TAKEN OUT OF ASTEC PACKAGE #######
###################################################################

#VIRTUAL EMBRYO
vo_login="emmanuel.faure"
vo_passwd="ascidie"
vo_lib_Path=astec_Path+"3DCloudEmbryo/" #Path for scripts library
vo_Path=postsegment_Path+"VO/" #Output path for mesh data
vo_files=vo_Path+EN+'_fuse_seg_post_vo_t$TIME.obj' # output files

