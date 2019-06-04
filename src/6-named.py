#!/usr/bin/python2.7


# Post correction of the segmentation
from definitions import *
if not os.path.isdir(postsegment_Path):
    os.mkdir(postsegment_Path) 
os.system("cp -f "+astec_Path+"definitions.py "+postsegment_Path )
os.system("cp -f "+astec_Path+"6-named.py "+postsegment_Path ) 


from lineage import write_tlp_from_lin_tree,read_lineage_tree,write_lineage_tree,timeNamed,timesNamed
from ascidians import process_multiple_contacts,name_propagation,fate_map
from lineage_test import pkl_lineage_test, imageDict

# MERGE Sophia - Montpellier (OPERATION A REFACTORISER DANS definitions.py)
# <<<<
# Names DEFINITION
Segmentation_method = raw_input("Select the nature of segmentation: 1 for ASTEC, 2 for GLASTEC\n")
if Segmentation_method =='1':
    print "ASTEC postcorrection..."
if Segmentation_method =='2':
    print "GLASTEC postcorrection..."
if (Segmentation_method =='2'):
	if not os.path.isdir(glastec_postsegment_Path):
		os.mkdir(glastec_postsegment_Path) 
	postsegment_files=glastec_postsegment_files
	post_lineage_tree_filename=glastec_post_lineage_tree_filename
# <<<<

lin_tree_information=read_lineage_tree(post_lineage_tree_filename) # Read the lineage tree (in case it was previously created)

#APPLY NAMED CELLS
print 'Read '+name_file
lin_tree_named={}
for l in open(name_file,'r'):
    tab=l.rstrip('\n').split(':')
    lin_tree_named[begin*10**4+int(tab[0])]=tab[1]

# MERGE Sophia - Montpellier
# <<<<
cell_named=name_propagation(lin_tree_information['lin_tree'], lin_tree_named, postsegment_files,begin, end, delta)
lin_tree_information['Names']=cell_named


"""
cell_named,errors=name_propagation(lin_tree_information['lin_tree'], lin_tree_named, postsegment_files,begin, end, delta)
lin_tree_information['Names']=cell_named
lin_tree_information['Error']=errors
"""
# <<<<


#CALCUL CELL SURFACE CONTACT
surfaces,barycenters,contact_surf=process_multiple_contacts(postsegment_files, begin, end,delta)
lin_tree_information['surface']=surfaces
lin_tree_information['cell_cell_surface_information']=contact_surf
lin_tree_information['barycenter']=barycenters


#CALCUL FATE MAP
fate,fate2=fate_map(lin_tree_information['lin_tree'],lin_tree_information['Names'],begin,end,delta)
lin_tree_information['fate']=fate
lin_tree_information['fate2']=fate2


### SAVE THE FINAL LINEAGE TREE
write_lineage_tree(named_lineage_tree_filename,lin_tree_information)
#Save Lineage in Tulip format (http://tulip.labri.fr/TulipDrupal/)
#write_tlp_from_lin_tree(named_lineage_tree_filename.replace('.pkl','.tlp'), lin_tree_information, ['volumes_information','cell_cell_surface_information','Names','fate','fate2'])  

### PROCESS LINEAGE TREE FILE VERIFICATION
print 'PROCESS LINEAGE TREE VERIFICATION'
image_dict_seg_post=imageDict(postsegment_files.replace("$TIME","*"))
rapport=pkl_lineage_test(lin_tree_information, image_dict_seg_post, file_out=named_lineage_tree_test_filename)

print 'LINEAGE TREE VERIFICATION DONE'

