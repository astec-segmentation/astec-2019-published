import sys,os

#Add the ASTEC Function
sys.path.append('ASTEC')

import numpy as np
from REGISTRATION import spatial_registration, associateRefFloLabels, symmetry_plane, sisters_association, propagate_all_correspondences, associateSegmentationsAtTimesFromLineageCorrespondences, readLUT, segmentationRelabellingAtTime, temporal_affine_registration

from morpheme_lineage import Morpheme_lineage
from morpheme_lineage import lineage_relabelled_point_cloud_at_time        # Pour l'ecriture des nuages de points relabellises a partir du lineage
from morpheme_lineage import labels_at_time_to_relabelled_correspondences  # Pour convertir les correspondances cell-snapshot to cell-snapshot etablies a l'aide de spatial_registration sur des labels originaux en correspondances de cellules relabellisees
from morpheme_lineage import write_correspondences  # Pour sauvegarder au format texte les dictionnaires de correspondances label-label (ou relabel-relabel)

#from morpheme_lineage import relabelled_to_labels_at_time_correspondences  # Complementaire de labels_at_time_to_relabelled_correspondences
root_data="/media/DATA/PROD/"
################
#### RALPH  ####
################

ralphKey='160707-Ralph-St8'

path_images_ralph_fused=root_data+"%s/FUSE/"%ralphKey
path_images_ralph_seg_post=root_data+"%s/FUSE/SEG/POST/"%ralphKey
path_lineage_ralph_seg_post=root_data+"%s/FUSE/SEG/POST/%s_fuse_seg_post_lin_tree.pkl"%(ralphKey,ralphKey)
path_lineage_ralph_seg_post_cor=root_data+"%s/FUSE/SEG/POST/%s_fuse_seg_post_lin_tree_corrected.pkl"%(ralphKey,ralphKey)
ralph_fused_ext='.inr'
ralph_seg_post_ext='.inr'

lineage_ralph_seg_post_cor = Morpheme_lineage( 'Ralph_seg_post', path_lineage_ralph_seg_post_cor ) 


#################
#### PATRICK ####
#################

patrickKey='140317-Patrick-St8'

path_images_patrick_fused=root_data+"%s/FUSE/"%patrickKey
path_images_patrick_seg_post=root_data+"%s/FUSE/SEG/POST/"%patrickKey
path_folder_patrick_reg=root_data+"%s/FUSE/REG/"%patrickKey
path_template_patrick_reg=root_data+"%s/FUSE/REG/%s_reg_compose_template.inr.gz"%(patrickKey,patrickKey)
path_lineage_patrick_seg_post_named=root_data+"%s/FUSE/SEG/POST/%s_fuse_seg_post_lin_tree_named.pkl"%(patrickKey,patrickKey)
path_lineage_patrick_seg_post_named_cor=root_data+"%s/FUSE/SEG/POST/%s_fuse_seg_post_lin_tree_named_corrected.pkl"%(patrickKey,patrickKey)
patrick_fused_ext='.inr.gz'
patrick_seg_post_ext='.inr.gz'

lineage_patrick_seg_post_named_cor = Morpheme_lineage( 'Patrick_seg_post_named', path_lineage_patrick_seg_post_named_cor )

###################
#### WORKSPACE ####
###################

WORKSPACE="WORKSPACE_example-registration/"
RESULTS="RESULTS_example-registration/"

##############################################################################################
########## spatial_registration applied at ref and flo timepoints of interest (toi) ##########
##############################################################################################
trsf_type='affine'

# Early registration (64 cells)
flo_time_point=4
ref_time_point=1

ralph_fused_toi_file = path_images_ralph_fused + ralphKey + '_fuse_t{0:03d}'.format(flo_time_point) + ralph_fused_ext
ralph_seg_post_toi_file = path_images_ralph_seg_post + ralphKey + "_fuse_seg_post_t{0:03d}".format(flo_time_point) + ralph_seg_post_ext
patrick_fused_toi_file = path_images_patrick_fused + patrickKey + '_fuse_t{0:03d}'.format(ref_time_point) + patrick_fused_ext
patrick_seg_post_toi_file = path_images_patrick_seg_post + patrickKey + "_fuse_seg_post_t{0:03d}".format(ref_time_point) + patrick_seg_post_ext

path_trsf_flo_ref=RESULTS+trsf_type+"_"+"Ralph_t{0:03d}_Patrick_t{1:03d}.trsf".format(flo_time_point, ref_time_point)
path_pairs_ref_flo=RESULTS+trsf_type+"_"+"Patrick_t{0:03d}_Ralph_t{1:03d}.pairs".format(ref_time_point, flo_time_point)

path_trsf_flo_ref, path_pairs_ref_flo, path_dices_ref_flo, path_residuals_ref_flo=spatial_registration(patrick_seg_post_toi_file, ralph_seg_post_toi_file, ref_fused_file=patrick_fused_toi_file, flo_fused_file=ralph_fused_toi_file,  # Input images
	path_trsf_flo_ref=path_trsf_flo_ref, path_pairs_ref_flo=path_pairs_ref_flo, # Output files
	folder_tmp=WORKSPACE, # Temporary files and paths
	trsf_type=trsf_type, estimator='lts', lts_fraction=0.9, # arguments for planeRegistration
	keep_mem=True, keep_bin=True, keep_hist=True, keep_sym=True, keep_inter=True, verbose=True)

path_patrick_toi_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_ref_associated.mha".format(ref_time_point, flo_time_point)
path_ralph_toi_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_flo_associated.mha".format(ref_time_point, flo_time_point)
path_labels_out_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_labels_associated.pairs".format(ref_time_point, flo_time_point)

path_patrick_toi_associated_cor=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_ref_associated_cor.mha".format(ref_time_point, flo_time_point)
path_ralph_toi_associated_cor=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_flo_associated_cor.mha".format(ref_time_point, flo_time_point)
associateRefFloLabels(patrick_seg_post_toi_file, ralph_seg_post_toi_file, path_pairs_ref_flo, path_patrick_toi_associated, path_ralph_toi_associated, 
	path_labels_out=path_labels_out_associated, path_trsf_flo_ref=path_trsf_flo_ref, ref=True, visu=True, verbose=True)
lut_early=readLUT(path_pairs_ref_flo)

# Lut manual corrections

path_labels_out_associated_cor=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_labels_associated.pairs.cor".format(ref_time_point, flo_time_point)
lut_cor_early=lut_early.copy()
# correction associations errors 
lut_cor_early[58]=13
lut_cor_early[83]=10
lut_cor_early[90]=18
# addition of missing associations
lut_cor_early[64]=11
lut_cor_early[67]=24
lut_cor_early[35]=62
lut_cor_early[30]=64
lut_cor_early[59]=35
lut_cor_early[45]=38
# saving corrected correspondences
write_correspondences(lut_cor_early, path_labels_out_associated_cor)
lut_cor_early=readLUT(path_labels_out_associated_cor)
associateRefFloLabels(patrick_seg_post_toi_file, ralph_seg_post_toi_file, path_labels_out_associated_cor, path_patrick_toi_associated_cor, path_ralph_toi_associated_cor, 
	path_labels_out=None, path_trsf_flo_ref=path_trsf_flo_ref, ref=True, visu=True, verbose=True)

#ref_deaths_cor_early=[]
#flo_deaths_cor_early=[]
#annotations_cor_early=[]
#for k, v in lut_cor_early.items():
#	ref_deaths_cor_early.append(lineage_patrick_seg_post_named.death(k, ref_time_point))
#	flo_deaths_cor_early.append(lineage_ralph_seg_post.death(v, flo_time_point))
#	annotations_cor_early.append("Patrick t{0:03d} l{1:03d}".format(ref_time_point, k) + " ; " + "Ralph t{0:03d} l{1:03d}".format(flo_time_point, v) )
#	print(str(k)+' '+str(lineage_patrick_seg_post_named.death(k, ref_time_point))+' ' + str(v)+' ' + str(lineage_ralph_seg_post.death(v, flo_time_point)))


if False:
	# Late registration (121 cells)
	ref_time_point=45
	flo_time_point=28

	ralph_fused_toi_file = path_images_ralph_fused + ralphKey + '_fuse_t{0:03d}'.format(flo_time_point) + ralph_fused_ext
	ralph_seg_post_toi_file = path_images_ralph_seg_post + ralphKey + "_fuse_seg_post_t{0:03d}".format(flo_time_point) + ralph_seg_post_ext
	patrick_fused_toi_file = path_images_patrick_fused + patrickKey + '_fuse_t{0:03d}'.format(ref_time_point) + patrick_fused_ext
	patrick_seg_post_toi_file = path_images_patrick_seg_post + patrickKey + "_fuse_seg_post_t{0:03d}".format(ref_time_point) + patrick_seg_post_ext

	path_trsf_flo_ref=RESULTS+trsf_type+'_'+"Ralph_t{0:03d}_Patrick_t{1:03d}.trsf".format(flo_time_point, ref_time_point)
	path_pairs_ref_flo=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}.pairs".format(ref_time_point, flo_time_point)

	path_trsf_flo_ref, path_pairs_ref_flo, path_dices_ref_flo, path_residuals_ref_flo=spatial_registration(patrick_seg_post_toi_file, ralph_seg_post_toi_file, ref_fused_file=patrick_fused_toi_file, flo_fused_file=ralph_fused_toi_file,  # Input images
		path_trsf_flo_ref=path_trsf_flo_ref, path_pairs_ref_flo=path_pairs_ref_flo, # Output files
		folder_tmp=WORKSPACE, # Temporary files and paths
		trsf_type=trsf_type, estimator='lts', lts_fraction=0.9, # arguments for planeRegistration
		keep_mem=True, keep_bin=True, keep_hist=True, keep_sym=True, keep_inter=True, verbose=True)

	path_patrick_toi_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_ref_associated.mha".format(ref_time_point, flo_time_point)
	path_ralph_toi_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_flo_associated.mha".format(ref_time_point, flo_time_point)
	path_labels_out_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_labels_associated.pairs".format(ref_time_point, flo_time_point)
	associateRefFloLabels(patrick_seg_post_toi_file, ralph_seg_post_toi_file, path_pairs_ref_flo, path_patrick_toi_associated, path_ralph_toi_associated, 
		path_labels_out=path_labels_out_associated, path_trsf_flo_ref=path_trsf_flo_ref, ref=True, visu=True, verbose=True)
	lut_late=readLUT(path_pairs_ref_flo)
	path_labels_out_associated_cor=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_labels_associated.pairs.cor".format(ref_time_point, flo_time_point)
	#lut_cor_late=readLUT(path_labels_out_associated_cor)

	#ref_deaths_late=[]
	#flo_deaths_late=[]
	#annotations_late=[]
	#for k, v in lut_late.items():
	#	ref_deaths_late.append(lineage_patrick_seg_post_named.death(k, ref_time_point))
	#	flo_deaths_late.append(lineage_ralph_seg_post.death(v, flo_time_point))
	#	annotations_late.append("Patrick t{0:03d} l{1:03d}".format(ref_time_point, k) + " ; " + "Ralph t{0:03d} l{1:03d}".format(flo_time_point, v) )
	#	print(str(k)+' '+str(lineage_patrick_seg_post_named.death(k, ref_time_point))+' ' + str(v)+' ' + str(lineage_ralph_seg_post.death(v, flo_time_point)))

	#ref_deaths_cor_late=[]
	#flo_deaths_cor_late=[]
	#annotations_cor_late=[]
	#for k, v in lut_cor_late.items():
	#	ref_deaths_cor_late.append(lineage_patrick_seg_post_named.death(k, ref_time_point))
	#	flo_deaths_cor_late.append(lineage_ralph_seg_post.death(v, flo_time_point))
	#	annotations_cor_late.append("Patrick t{0:03d} l{1:03d}".format(ref_time_point, k) + " ; " + "Ralph t{0:03d} l{1:03d}".format(flo_time_point, v) )
	#	print(str(k)+' '+str(lineage_patrick_seg_post_named.death(k, ref_time_point))+' ' + str(v)+' ' + str(lineage_ralph_seg_post.death(v, flo_time_point)))




	# End gastrulation registration (182 cells)
	ref_time_point=54
	flo_time_point=38

	ralph_fused_toi_file = path_images_ralph_fused + ralphKey + '_fuse_t{0:03d}'.format(flo_time_point) + ralph_fused_ext
	ralph_seg_post_toi_file = path_images_ralph_seg_post + ralphKey + "_fuse_seg_post_t{0:03d}".format(flo_time_point) + ralph_seg_post_ext
	patrick_fused_toi_file = path_images_patrick_fused + patrickKey + '_fuse_t{0:03d}'.format(ref_time_point) + patrick_fused_ext
	patrick_seg_post_toi_file = path_images_patrick_seg_post + patrickKey + "_fuse_seg_post_t{0:03d}".format(ref_time_point) + patrick_seg_post_ext

	path_trsf_flo_ref=RESULTS+trsf_type+'_'+"Ralph_t{0:03d}_Patrick_t{1:03d}.trsf".format(flo_time_point, ref_time_point)
	path_pairs_ref_flo=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}.pairs".format(ref_time_point, flo_time_point)

	path_trsf_flo_ref, path_pairs_ref_flo, path_dices_ref_flo, path_residuals_ref_flo=spatial_registration(patrick_seg_post_toi_file, ralph_seg_post_toi_file, ref_fused_file=patrick_fused_toi_file, flo_fused_file=ralph_fused_toi_file,  # Input images
		path_trsf_flo_ref=path_trsf_flo_ref, path_pairs_ref_flo=path_pairs_ref_flo, # Output files
		folder_tmp=WORKSPACE, # Temporary files and paths
		trsf_type=trsf_type, estimator='lts', lts_fraction=0.9, # arguments for planeRegistration
		keep_mem=True, keep_bin=True, keep_hist=True, keep_sym=True, keep_inter=True, verbose=True)

	path_patrick_toi_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_ref_associated.mha".format(ref_time_point, flo_time_point)
	path_ralph_toi_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_flo_associated.mha".format(ref_time_point, flo_time_point)
	path_labels_out_associated=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_labels_associated.pairs".format(ref_time_point, flo_time_point)
	associateRefFloLabels(patrick_seg_post_toi_file, ralph_seg_post_toi_file, path_pairs_ref_flo, path_patrick_toi_associated, path_ralph_toi_associated, 
		path_labels_out=path_labels_out_associated, path_trsf_flo_ref=path_trsf_flo_ref, ref=True, visu=True, verbose=True)
	lut_gastrul=readLUT(path_pairs_ref_flo)
	#path_labels_out_associated_cor=RESULTS+trsf_type+'_'+"Patrick_t{0:03d}_Ralph_t{1:03d}_labels_associated.pairs.cor".format(ref_time_point, flo_time_point)
	#lut_cor_gastrul=readLUT(path_labels_out_associated_cor)

	#ref_deaths_gastrul=[]
	#flo_deaths_gastrul=[]
	#annotations_gastrul=[]
	#for k, v in lut_gastrul.items():
	#	ref_deaths_gastrul.append(lineage_patrick_seg_post_named.death(k, ref_time_point))
	#	flo_deaths_gastrul.append(lineage_ralph_seg_post.death(v, flo_time_point))
	#	annotations_gastrul.append("Patrick t{0:03d} l{1:03d}".format(ref_time_point, k) + " ; " + "Ralph t{0:03d} l{1:03d}".format(flo_time_point, v) )
	#	print(str(k)+' '+str(lineage_patrick_seg_post_named.death(k, ref_time_point))+' ' + str(v)+' ' + str(lineage_ralph_seg_post.death(v, flo_time_point)))

	#ref_deaths_cor_gastrul=[]
	#flo_deaths_cor_gastrul=[]
	#annotations_cor_gastrul=[]
	#for k, v in lut_cor_gastrul.items():
	#	ref_deaths_cor_gastrul.append(lineage_patrick_seg_post_named.death(k, ref_time_point))
	#	flo_deaths_cor_gastrul.append(lineage_ralph_seg_post.death(v, flo_time_point))
	#	annotations_cor_gastrul.append("Patrick t{0:03d} l{1:03d}".format(ref_time_point, k) + " ; " + "Ralph t{0:03d} l{1:03d}".format(flo_time_point, v) )
	#	print(str(k)+' '+str(lineage_patrick_seg_post_named.death(k, ref_time_point))+' ' + str(v)+' ' + str(lineage_ralph_seg_post.death(v, flo_time_point)))

	#for k,v in lut_early.items():
	#	if lineage_patrick_seg_post_named.death(k, ref_time_point) <

###########################################################################################################################
####### Correction des fichiers de lignage de Patrick et Ralph par rapport aux donnees de post-segmentation reelles #######
###########################################################################################################################

if False: 

	from lineage import read_lineage_tree
	from pkl_converter import getLongIDFromTimeAndLabel
	lin_ralph=read_lineage_tree(path_lineage_ralph_seg_post)
	lin_pat=read_lineage_tree(path_lineage_patrick_seg_post_named)


	bary_factor=0.3
	vol_factor=bary_factor**3

	lin_ralph['voxel_volume']={}
	lin_ralph['real_volume']={}
	lin_ralph['voxel_barycenter']={}
	lin_ralph['real_barycenter']={}

	for time in range(4, 96):
		cloud_file=path_images_ralph_seg_post+'CLOUD/'+ ralphKey + "_fuse_seg_post_t{0:03d}".format(time) + ".cloud"
		cloud=np.loadtxt(cloud_file,comments="#",unpack=False, delimiter=" ",skiprows=1)
		size=cloud.shape

		for i in range(size[0]):

			lin_ralph['real_barycenter'][getLongIDFromTimeAndLabel(time, int(cloud[i,0]))] = (bary_factor*cloud[i,1],bary_factor*cloud[i,2],bary_factor*cloud[i,3])
			lin_ralph['real_volume'][getLongIDFromTimeAndLabel(time, int(cloud[i,0]))] = vol_factor*cloud[i,4]
			lin_ralph['voxel_barycenter'][getLongIDFromTimeAndLabel(time, int(cloud[i,0]))] = (cloud[i,1],cloud[i,2],cloud[i,3])
			lin_ralph['voxel_volume'][getLongIDFromTimeAndLabel(time, int(cloud[i,0]))] = cloud[i,4]

	lin_pat[1]['voxel_volume']={}
	lin_pat[1]['real_volume']={}
	lin_pat[1]['voxel_barycenter']={}
	lin_pat[1]['real_barycenter']={}

	for time in range(1, 193):
		cloud_file=path_images_patrick_seg_post+'CLOUD/'+ patrickKey + "_fuse_seg_post_t{0:03d}".format(time) + ".cloud"
		cloud=np.loadtxt(cloud_file,comments="#",unpack=False, delimiter=" ",skiprows=1)
		size=cloud.shape

		for i in range(size[0]):

			lin_pat[1]['real_barycenter'][getLongIDFromTimeAndLabel(time, int(cloud[i,0]))] = (bary_factor*cloud[i,1],bary_factor*cloud[i,2],bary_factor*cloud[i,3])
			lin_pat[1]['real_volume'][getLongIDFromTimeAndLabel(time, int(cloud[i,0]))] = vol_factor*cloud[i,4]
			lin_pat[1]['voxel_barycenter'][getLongIDFromTimeAndLabel(time, int(cloud[i,0]))] = (cloud[i,1],cloud[i,2],cloud[i,3])
			lin_pat[1]['voxel_volume'][getLongIDFromTimeAndLabel(time, int(cloud[i,0]))] = cloud[i,4]

	import cPickle as pickle

	output = open(path_lineage_ralph_seg_post_cor, 'wb')
	pickle.dump(lin_ralph, output)
	output.close()

	output = open(path_lineage_patrick_seg_post_named_cor, 'wb')
	pickle.dump(lin_pat, output)
	output.close()

####################################################
#### GENERATE POINTS CLOUDS WITH RELABELLED IDS ####
####################################################

if False:
	for time in range(4,lineage_ralph_seg_post_cor.timepoints()[-1]+1):
		lineage_relabelled_point_cloud_at_time(lineage_ralph_seg_post_cor, time, '/home/gmicheli/DIG-EM/Data/160707-Ralph-St8/FUSE/SEG/POST/CLOUD/160707-Ralph-St8_fuse_seg_post_relabelled_t%03d.cloud'%time)

	for time in range(1, lineage_patrick_seg_post_named_cor.timepoints()[-1]+1):
		lineage_relabelled_point_cloud_at_time(lineage_patrick_seg_post_named_cor, time, '/home/gmicheli/DIG-EM/Data/140317-Patrick-St8/FUSE/SEG/POST/CLOUD/140317-Patrick-St8_fuse_seg_post_relabelled_t%03d.cloud'%time)

########################
#### CORRESPONDENCES ###
########################

from cpp_wrapping import compose_trsf

#### CORRESPONDENCES READING AND PROPAGATION
time_point_ref=1
time_point_flo=4

#path_template_patrick_reg="/home/gmicheli/DIG-EM/Data/140317-Patrick-St8/FUSE/REG/140317-Patrick-St8_reg_compose_template.inr.gz"

#file_correspondences='affine_Patrick_t%03d_Ralph_t%03d_labels_associated.pairs.cor'%(time_point_ref, time_point_flo)

#lut_cor_early=readLUT(path_labels_out_associated_cor)
# On convertit les correspondences etablies a partir des labels originaux des images segmentees a des instants specifiques en correspondences de cellules relabellisees
correspondences_init=labels_at_time_to_relabelled_correspondences(lut_cor_early,
															 lineage_patrick_seg_post_named_cor, lineage_ralph_seg_post_cor, time_point_ref, time_point_flo)

final_correspondences={}
unpropagated_ref={}
scalars={}
file_new_correspondences={}
a={}
b={}
for threshold_angle in [45, 60, 90, 30, 40]:

	file_new_correspondences[threshold_angle]=RESULTS+"correspondences_from_early_cor_%d_degrees_relabel.lut"%threshold_angle
	final_correspondences[threshold_angle], unpropagated_ref[threshold_angle], scalars[threshold_angle]=propagate_all_correspondences(lineage_patrick_seg_post_named_cor, lineage_ralph_seg_post_cor, correspondences_init, threshold_angle=threshold_angle, verbose=True )
	write_correspondences(final_correspondences[threshold_angle], file_new_correspondences[threshold_angle], ID_key='RELABEL')

	if True:
		for time_ref in lineage_patrick_seg_post_named_cor.timepoints():
			path_seg_post_patrick_at_time = path_images_patrick_seg_post+patrickKey+'_fuse_seg_post_t%03d'%time_ref+patrick_seg_post_ext
			path_trsf_at_time = path_folder_patrick_reg+patrickKey+"_reg_compose_t%03d.trsf"%time_ref
			path_ref_out=RESULTS+"RELABELLED/relabelled_Patrick_t%03d_aligned_with_Ralph_angle_%d.mha"%(time_ref,threshold_angle)

			segmentationRelabellingAtTime(path_seg_post_patrick_at_time, path_ref_out, lineage_patrick_seg_post_named_cor, time_ref, lut=final_correspondences[threshold_angle].keys(), trsf=path_trsf_at_time, template=path_template_patrick_reg,
										 voxelsize=None, iso=None, dimensions=None, backgroundLabels=[0,1], visu=True, verbose=True)
		for time_flo in lineage_ralph_seg_post_cor.timepoints():
			path_seg_post_ralph_at_time = path_images_ralph_seg_post+ralphKey+'_fuse_seg_post_t%03d.'%time_flo+ralph_seg_post_ext
			a[threshold_angle],b[threshold_angle]=temporal_affine_registration(lineage_patrick_seg_post_named_cor, lineage_ralph_seg_post_cor, final_correspondences[threshold_angle], include_end_of_time=False)
			time_ref=int((time_flo-b[threshold_angle])/a[threshold_angle]) # temporal registration given by linear regression
			if time_ref<lineage_patrick_seg_post_named_cor.timepoints()[0]:
				time_ref=lineage_patrick_seg_post_named_cor.timepoints()[0]
			if time_ref>lineage_patrick_seg_post_named_cor.timepoints()[-1]:
				lineage_patrick_seg_post_named_cor.timepoints()[-1]
			print "TIME_REF = t%03d, TIME_FLO = t%03d"%(time_ref,time_flo)
			path_flo_out=RESULTS+"RELABELLED/relabelled_Ralph_t%03d_on_Patrick_t%03d_angle_%d.mha"%(time_flo,time_ref, threshold_angle)
			# Calcul de la transformation
			out=pointCloudRegistration(lineage_patrick_seg_post_named_cor.relabelled_barycenters_at_time(time_ref), lineage_ralph_seg_post_cor.relabelled_barycenters_at_time(time_flo), 
					final_correspondences[threshold_angle], skip_not_found=True, lazy=False, verbose=True)
			trsf_f_r=out['trsf']
			#path_trsf_r_template = "/home/gmicheli/DIG-EM/Data/140317-Patrick-St8/FUSE/REG/140317-Patrick-St8_reg_compose_t%03d.trsf"%time_ref
			path_trsf_at_time = path_folder_patrick_reg+patrickKey+"_reg_compose_t%03d.trsf"%time_ref
			path_trsf_tmp = "tmp_path.trsf"
			f=open(path_trsf_tmp,"w")
			f.write(str(np.array(trsf_f_r)).replace('[','').replace(']',''))
			f.close()
			# Composition de la transformation avec la transformation de la reference vers le template
			compose_trsf(path_trsf_tmp, path_trsf_at_time, path_trsf_tmp)
	
			segmentationRelabellingAtTime(path_seg_post_ralph_at_time, path_flo_out, lineage_ralph_seg_post_cor, time_flo, reverse_lut=final_correspondences[threshold_angle], trsf=path_trsf_tmp, template=path_template_patrick_reg,
					voxelsize=None, iso=None, dimensions=None, backgroundLabels=[0,1], visu=True, verbose=True)

			if os.path.exists(path_trsf_tmp):
				os.system('rm %s'%path_trsf_tmp)

