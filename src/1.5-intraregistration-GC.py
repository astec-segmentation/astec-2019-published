from definitions import *
import numpy as np
import scipy
from scipy import *
import cPickle as pkl
import os
import sys
from ImageHandling import imread, imsave, SpatialImage
from lineage import read_lineage_tree,timeNamed
from computebarycenters import Thread_BaryCenter, WaitThreadsEnd

print "/!\\ Caution, this script has not been officially integrated yet and should be reviewed /!\\ "

if not os.path.isdir(reg_glas_Path):
    os.mkdir(reg_glas_Path)  
os.system("cp -f "+astec_Path+"definitions.py "+reg_glas_Path )
os.system("cp -f "+astec_Path+"7-registration.py "+reg_glas_Path )

def save_pkl(obj, filename):
    #saves pickle object, with name given by "filename".
    with open(filename, 'wb') as output:
        pkl.dump(obj, output, pkl.HIGHEST_PROTOCOL)

def load_pkl(filename):
    #loads pickle object
    with open(filename, 'r') as file:
        obj=pkl.load(file)
    return(obj)

#Time to register:
t_in=begin
end=95
t_fin=end


def save_object(obj, filename):
	#saves pickle object, with name given by "filename".
    with open(filename, 'wb') as output:
        pkl.dump(obj, output, pkl.HIGHEST_PROTOCOL)

def load_object(filename):
	#loads pickle object
    with open(filename, 'r') as file:
        obj=pkl.load(file)
   	return(obj)

def diags(i, Q):
	return (i**2)*(1-cos(Q))+cos(Q)

def off_diags(i, j, k, p, Q):
	return (1-cos(Q))*i*j+p*k*sin(Q)

def unit_vec(v):
	n=sqrt(np.dot(v,v))
	return [(1.0/n)*k for k in v]

def proj_perp(vec, perp):
	v1=(np.dot(np.array(vec),np.array(perp)))*np.array(unit_vec(perp))
	v2=np.array(vec)-v1
	if np.dot(v2,vec)<0:
		v2=-1.0*np.array(v2)
	return unit_vec(v2)

def rot_mat(Q, ax):
	axx=unit_vec(ax)
	i=axx[0]
	j=axx[1]
	k=axx[2]
	return [[diags(i, Q), off_diags(i,j,k,-1.0,Q), off_diags(i,k,j,+1.0,Q)],[off_diags(j,i,k,+1.0,Q), diags(j,Q), off_diags(j,k,i,-1.0,Q)],[off_diags(k,i,j,-1.0,Q),off_diags(k,j,i,+1.0,Q),diags(k,Q)]]

def Rotation_Matrices(angles, axes, t_start):
	barys_in=[bary[k] for k in time_labels[t_start]]
	rot_center=np.mean(barys_in,axis=0)
	T_back=-1.0*np.array(rot_center)
	Rot=[rot_mat(angles[k], axes[k]) for k in range(0,len(angles))]
	M=Rot[0]
	for p in range(1,len(Rot)):
		M=np.dot(Rot[p],M)
	T=np.dot(M,T_back)-T_back
	return M, T

def saveMatrix(file, filename):
	while os.path.exists(filename)==True:
		ans=" "
		while ans not in ["y","n"]:
			ans=raw_input("File already exists! Do you want to overwrite? (y/n) ")
		if ans=="y":
			os.remove(filename)
			break
		else:
			filename=raw_input("Please enter a new file name: ")
	if "." not in filename:
		filename=filename+".txt"
	with open(filename, "a") as f:
		for i in file:
			f.write(str(i[0])+"\t"+str(i[1])+"\t"+str(i[2])+"\t"+str(i[3])+"\n")
	print "File saved succesfully!"
	print ""
	return None



def saveDictPoints(file, filename, selection, size):
	while os.path.exists(filename)==True:
		ans=" "
		while ans not in ["y","n"]:
			ans=raw_input("File already exists! Do you want to overwrite? (y/n) ")
		if ans=="y":
			os.remove(filename)
			break
		else:
			filename=raw_input("Please enter a new file name: ")
	if "." not in filename:
		filename=filename+".txt"
	with open(filename, "a") as f:
		for i in file.keys():
			f.write(str(i)+":"+str(file[i][0])+","+str(file[i][1])+","+str(file[i][2])+":"+str(size)+":"+str(selection)+"\n")
	print "File saved succesfully!"
	print ""
	return None


def saveDictPointsVols(file, filename, selection, size_ratio, ref_vol=60000):
	while os.path.exists(filename)==True:
		ans=" "
		while ans not in ["y","n"]:
			ans=raw_input("File already exists! Do you want to overwrite? (y/n) ")
		if ans=="y":
			os.remove(filename)
			break
		else:
			filename=raw_input("Please enter a new file name: ")
			if "." not in filename:
				filename=filename+".txt"
	with open(filename, "a") as f:
		for i in file.keys():
			size=size_ratio+0.03*size_ratio*1.0*vol[i]/ref_vol
			f.write(str(i)+":"+str(file[i][0])+","+str(file[i][1])+","+str(file[i][2])+":"+str(size)+":"+str(selection)+"\n")
	print "File saved succesfully!"
	print ""
	return None

def point_sets(t_start):
	lab_in=time_labels[t_start]
	lab_out=time_labels[t_start+1]
	set_in=[k for k in lab_in if len(lin.get(k,[3,3,2]))<1.5 and type(bary.get(k,"NO")) is not str and type(bary.get(lin[k][0],"NO")) is not str]
	p_t=[bary[k] for k in set_in]
	p_tp1=[bary[lin[k][0]] for k in set_in]
	return p_t,p_tp1

def rigid_transform_3D(A, B):
	assert len(A) == len(B)
	N = len(A) # total points
	centroid_A = mean(A, axis=0)
	centroid_B = mean(B, axis=0)
	H = np.sum([np.outer(np.array(A[k])-np.array(centroid_A),np.array(B[k])-np.array(centroid_B)) for k in range(0,len(A))],axis=0)
	U, S, Vt = np.linalg.svd(H)
	R = np.dot(transpose(Vt),transpose(U))
	if np.linalg.det(R) < 0:
		print "Reflection detected"
		Vt[2,:] *= -1
		R =  np.dot(transpose(Vt),transpose(U))
	t = -np.dot(R,centroid_A) + centroid_B
	return R, t

def saveAxes(file, filename):
	while os.path.exists(filename)==True:
		ans=" "
		while ans not in ["y","n"]:
			ans=raw_input("File already exists! Do you want to overwrite? (y/n) ")
		if ans=="y":
			os.remove(filename)
			break
		else:
			filename=raw_input("Please enter a new file name: ")
			if "." not in filename:
				filename=filename+".txt"
	with open(filename, "a") as f:
		f.write(str(20033)+":"+str(bary_in[0])+","+str(bary_in[1])+","+str(bary_in[2])+":"+str(x_ax[0])+","+str(x_ax[1])+","+str(x_ax[2])+":0.3"+":0.1"+":1"+"\n")
		f.write(str(20033)+":"+str(bary_in[0])+","+str(bary_in[1])+","+str(bary_in[2])+":"+str(y_ax[0])+","+str(y_ax[1])+","+str(y_ax[2])+":0.3"+":0.1"+":3"+"\n")
		f.write(str(20033)+":"+str(bary_in[0])+","+str(bary_in[1])+","+str(bary_in[2])+":"+str(z_ax[0])+","+str(z_ax[1])+","+str(z_ax[2])+":0.3"+":0.1"+":7"+"\n")
	print "File saved succesfully!"
	print ""
	return None

def full_transf_matrix(Rot, Trasl):
	FM=[]
	for i in range(0,3):
		FMi=[]
		for r in Rot[i]:
			FMi.append(r)
		FMi.append(Trasl[i])
		FM.append(FMi)
	FM.append([0,0,0,1])
	return FM


def comparison(Rot, Trasl, point):
	FM=full_transf_matrix(Rot, Trasl)
	p1=np.dot(Rot, point)+Trasl
	p2a=[k for k in point]
	p2a.append(1)
	p2=np.dot(FM, p2a)
	return p1, p2


def prod_list_matrices(mats):
	Mt=mats[0]
	for x in range(1,len(mats)):
		Mt=np.dot(Mt,mats[x])
	return Mt

def compose_transformation_back(tf, t0):
	times=range(t0,tf+1)
	mats=[Mm[t] for t in times]
	Mcomp=prod_list_matrices(mats)
	partials=np.sum([np.dot(prod_list_matrices(mats[0:i]),Tm[times[i]]) for i in range(1,len(mats))],axis=0)
	Tcomp=partials+np.array(Tm[t0])
	return Mcomp, Tcomp

def transform_point_compos(point, Mc, Tc):
	pt=np.dot(Mc,point)+np.array(Tc)
	return [int(round(k)) for k in pt]


def proj_perp(vec, perp):
	v1=(np.dot(np.array(vec),np.array(perp)))*np.array(unit_vec(perp))
	v2=np.array(vec)-v1
	if np.dot(v2,vec)<0:
		v2=-1.0*np.array(v2)
	return unit_vec(v2)











#Calcul barycenters
"""
if not os.path.isfile(glastec_postsegment_Path+'new_barycenters.pkl'):
    print 'Compute Barycenters for ' + EN
    Threads=[]
    for t in range(begin,end+1):
        if not os.path.isfile(glastec_postsegment_Path+'new_barycenters_t'+str(t)+'.pkl'):
            th=Thread_BaryCenter(timeNamed(glastec_postsegment_files,t),t)
            th.start()
            Threads.append(th)
        
        if len(Threads)>35:
            WaitThreadsEnd(Threads)
            Threads=[]

    #Merge...
    bary={}
    for t in range(begin,end+1):
            baryT=load_pkl(glastec_postsegment_Path+'new_barycenters_t'+str(t)+'.pkl')
            print 'Found '+str(len(baryT))+' at '+str(t)
            for c in baryT:
                bary[c]=baryT[c]

    save_pkl(bary,glastec_postsegment_Path+'new_barycenters.pkl')

"""

lin_tree_information=read_lineage_tree(glastec_post_lineage_tree_filename) # Read the lineage tree (in case it was previously created)
#barycenter
#bary=lin_tree_information['barycenter']
bary=load_pkl(glastec_postsegment_Path+'new_barycenters_1.pkl')



#lineage
print lin_tree_information.keys()
lin=lin_tree_information['lin_tree']
volume= lin_tree_information['volumes_information']

#print volume[10083]
#names
print 'Read '+name_file
lin_tree_named={}
for l in open(name_file,'r'):
    tab=l.rstrip('\n').split(':')
    lin_tree_named[begin*10**4+int(tab[0])]=tab[1]
print lin_tree_named

names=lin_tree_named


#names=load_object("names")
#ROTATIONS NEEDED TO BRING GERM CELLS DOWN AT THE POINT t_in: SPECIFY AXES OF ROTATION AND ANGLES. TO ADJUST MANUALLY FOR EACH DATASET#

cells=bary.keys()
time_labels={}
for t in range(t_in,t_fin+1):
	labs_t=[k for k in cells if k/10**4==t]
	time_labels.update({t:labs_t})


print [k for k in bary.keys() if k not in names.keys() and k in time_labels[1]]
print [k for k in names.keys() if k not in bary.keys()]


cgerms=[k for k in time_labels[t_in] if "b7.0006" in names[k]]
ctop=[k for k in time_labels[t_in] if "a7.0010" in names[k]]

cveg=[k for k in time_labels[t_in] if "a7.0002" in names[k]]
can=[k for k in time_labels[t_in] if "a7.0012" in names[k] or "a7.0016" in names[k]]

btop=np.mean([bary[k] for k in ctop],axis=0)
bgerms=np.mean([bary[k] for k in cgerms],axis=0)

bveg=np.mean([bary[k] for k in cveg],axis=0)
ban=np.mean([bary[k] for k in can],axis=0)

dir_tg=unit_vec(btop-bgerms)
goal_dir_1=[1,0,0]
ang_tg_x=arccos(np.dot(goal_dir_1,dir_tg))
rot_ax_1=unit_vec(np.cross(dir_tg,goal_dir_1))

mm=rot_mat(ang_tg_x, rot_ax_1)

dir_va=np.dot(mm,unit_vec(bveg-ban))
dir_va=proj_perp(dir_va, goal_dir_1)
goal_dir_2=[0,0,-1]
ang_va_z=arccos(np.dot(goal_dir_2,dir_va))
rot_ax_2=unit_vec(np.cross(dir_va, goal_dir_2))

angs=[ang_tg_x,ang_va_z]

axs=[rot_ax_1,rot_ax_2]

M_g,T_g=Rotation_Matrices(angs,axs,t_in) #INITIAL ROTATIONS FOR GERM CELLS DOWN

Mp={}
Mm={}
Tp={}
Tm={}
for t in range(t_in,t_fin):
	pt,ptp1=point_sets(t)
	Mf,Tf=rigid_transform_3D( pt, ptp1 )
	Mb,Tb=rigid_transform_3D( ptp1, pt )
	Mp.update({t:Mf})
	Mm.update({t+1:Mb})
	Tp.update({t:Tf})
	Tm.update({t+1:Tb})
Mm.update({t_in:M_g})
Tm.update({t_in:T_g})

save_object(Mm,"Mback")
save_object(Mp,"Mforw")


















# Calcul registration files (Bruno)

comp_Mback={}
comp_Tback={}
for t in range(t_in,t_fin+1):
	Mc,Tc=compose_transformation_back(t, t_in)
	comp_Mback.update({t:Mc})
	comp_Tback.update({t:Tc})


save_object(comp_Mback,reg_glas_Path+EN+'_Mback.pkl')
save_object(comp_Tback,reg_glas_Path+EN+'_Tback.pkl')


bary_reg={}
for t in range(t_in,t_fin+1):
	cells=time_labels[t]
	Mc,Tc=compose_transformation_back(t, t_in)
	bary_reg.update({k:np.array(transform_point_compos(bary[k], Mc, Tc)) for k in cells if type(bary.get(k,"NO")) is not str})










# Apply the transformation to post segmentation files (Emmanuel)

print 'Load Transformation '
comp_Mback=load_pkl(reg_glas_Path+EN+'_Mback.pkl')
comp_Tback=load_pkl(reg_glas_Path+EN+'_Tback.pkl')

identity=np.zeros([3,3])
for i in range(3):
    identity[i,i]=1.0


#First We compte the largest Transformed Images
minB=[100000,100000,100000]
maxB=[-100000,-100000,-100000]
if not os.path.isfile(reg_glas_Path+EN+'_MinMaxTransformed.pkl') :
    print 'Compte Max Image Size'
    for t in range(begin,end+1):
        print 'Read '+timeNamed(glastec_postsegment_files,t)
        segmented_img  = imread(timeNamed(glastec_postsegment_files,t))
       
        Mback=identity
        Tback=np.zeros(3)
        if t in comp_Mback and t in comp_Tback:
            Mback=comp_Mback[t]
            Tback=comp_Tback[t]

        def trans(coords):
            t_coords=np.dot(Mback,coords) #Rotation
            new_coords=np.zeros_like(coords)
            for i in range(3):
                new_coords[i]=np.round(t_coords[i]+Tback[i]).astype(np.int64) #Translation
            
            for i in range(3):
                minB[i]=min(minB[i],new_coords[i].min())
                maxB[i]=max(maxB[i],new_coords[i].max())
            #print "minB="+str(minB)+" maxB="+str(maxB)


        Xmax=segmented_img.shape[0]
        Ymax=segmented_img.shape[1]
        Zmax=segmented_img.shape[2]

        trans([0,0,0])
        trans([0,Ymax,0])
        trans([0,Ymax,Zmax])
        trans([0,0,Zmax])
        trans([Xmax,0,0])
        trans([Xmax,Ymax,0])
        trans([Xmax,Ymax,Zmax])
        trans([Xmax,0,Zmax])

        print "minB="+str(minB)+" maxB="+str(maxB)
    save_pkl([minB,maxB],reg_glas_Path+EN+'_MinMaxTransformed.pkl')
else:
    [minB,maxB]=load_pkl(reg_glas_Path+EN+'_MinMaxTransformed.pkl')
    print "minB="+str(minB)+" maxB="+str(maxB)


#The new Image Shape
new_shape=np.zeros(3, dtype=np.int)
for i in range(3):
    new_shape[i]=int(round(1+maxB[i]-minB[i]))
print 'New Shape '+str(new_shape)

#postsegment_files="/Users/emmanuelfaure/Projects/160606-4DCloudEmbryos/ConvertEmbryoToMesh/170511-Emmanuel-St8/SEG/"+EN+"_fuse_seg_post_t$TIME.inr"
#reg_files="/Users/emmanuelfaure/Projects/160606-4DCloudEmbryos/ConvertEmbryoToMesh/170511-Emmanuel-St8/"+EN+"_fuse_seg_post_reg_t$TIME.inr"
#We apply this transformation
for t in range(begin,end+1):
    reg_file=timeNamed(reg_glas_files,t)
    print 'Process '+reg_file
    if not os.path.isfile(reg_file):
        print 'Read '+timeNamed(glastec_postsegment_files,t)
        segmented_img  = imread(timeNamed(glastec_postsegment_files,t))

        Mback=identity
        Tback=np.zeros(3)
        if t in comp_Mback and t in comp_Tback:
            print 'Found transformation at '+str(t)
            Mback=comp_Mback[t]
            Tback=comp_Tback[t]

        coords=np.where(segmented_img>1) #1 IS BACKGROUND
        pixels=segmented_img[coords]

        t_coords=np.dot(Mback,coords) #Rotation
        new_coords=np.zeros_like(coords)
        for i in range(3):
            new_coords[i]=np.round(t_coords[i]+Tback[i]-minB[i]).astype(np.int64) #Translation

        new_img=SpatialImage(np.zeros([new_shape[0],new_shape[1],new_shape[2]],dtype=np.uint16))

        new_img[[new_coords[0],new_coords[1],new_coords[2]]]=pixels
        

        img_to_dilate=np.zeros([new_shape[0],new_shape[1],new_shape[2]],dtype=np.uint16)
        img_to_dilate[[new_coords[0],new_coords[1],new_coords[2]]]=1
        img_dilate=scipy.ndimage.binary_dilation(img_to_dilate).astype(img_to_dilate.dtype)
        img_dilate=img_dilate-img_to_dilate
        del img_to_dilate
        #coords=np.where(img_dilate==1)
        #del img_dilate


        #pixels=segmented_img[prev_coords]
        #new_img[coords]=pixels#segmented_img[prev_coords]

        #del img_dilate
        
        for z in range(new_shape[2]):
            coords=np.where(img_dilate[:,:,z]==1)
            for ixy in range(len(coords[0])):
                x=coords[0][ixy]
                y=coords[1][ixy]
                allfill=True
                pixelsValue=0
                for xx in range(x-1,x+1):
                    if allfill:
                        for yy in range(y-1,y+1):
                            if allfill and not (xx==x and yy==y):
                                if new_img[xx,yy,z]==0:
                                    allfill=False
                                else:
                                    pixelsValue=new_img[xx,yy,z]
                    if allfill:
                        new_img[x,y,z]=pixelsValue

                
        imsave(reg_file,new_img)
