import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__),"CommunFunctions"))
from ImageHandling import imread, imsave, SpatialImage
from scipy import ndimage as nd

from ACE import GLACE_from_resampled_segmentation, GACE
from cpp_wrapping import (obsolete_non_linear_registration, apply_trsf, obsolete_watershed,
                          obsolete_find_local_minima, outer_detection,
                          obsolete_gradient_norm, obsolete_morpho, obsolete_copy, obsolete_mc_adhocFuse, obsolete_Arit)


def compute_volumes(im, labels = None, real = True):
    """
    Return a dictionary, { label: volume }
    im : SpatialImage (label image)
    labels : list of labels for which to compute the volumes (if None, all volumes are computed)
    """
    labels = np.unique(im)

    volume = nd.sum(np.ones_like(im), im, index=np.int16(labels))
    return dict(zip(labels, volume))


def create_seed(parameters): 
    """
    Erodes the label i in the binary image tmp
    tmp : binary SpatialImage
    max_size_cell : size max allow for a cell (here put at np.inf)
    size_cell : size of the cell to erode
    iterations : maximum number of iterations for normal cells
    out_iterations : maximum number of iterations for exterior
    bb : bounding box if tmp in the global image (necessary when computing in parallel)
    i : label of the cell to erode
    """
    tmp, max_size_cell, size_cell, iterations, out_iterations, bb, i=parameters
    nb_iter=iterations
    if i==1:
        nb_iter=out_iterations
        opened=nd.binary_erosion(tmp, iterations=nb_iter)
        while len(nd.find_objects(opened))!=1 and nb_iter>=0:
            nb_iter-=1
            opened=nd.binary_erosion(tmp, iterations=nb_iter)
    else:
        opened=nd.binary_erosion(tmp, iterations=nb_iter)
        while len(nd.find_objects(opened))!=1 and nb_iter>=0:
            nb_iter-=1
            opened=nd.binary_erosion(tmp, iterations=nb_iter)
    if max_size_cell<size_cell:
        num=1
    else:
        num=i
    return opened, num, bb
    

def create_seeds(seg, max_size_cell=np.inf, min_size_cell=1000, iterations=10, out_iterations=25, nb_proc=26):
    """
    Erodes all the labels in the segmented image seg
    seg : Segmentation to erode (SpatialImage)
    max_size_cell : size maximum of a cell in number of voxels
    min_size_cell : size minimum of a cell in number of voxels
    iterations : maximum number of iterations for normal cells
    out_iterations : maximum number of iterations for exterior
    nb_proc : number maximum of processors allowed to be used
    """
    from multiprocessing import Process, Queue, Pool
    bboxes=nd.find_objects(seg)
    seeds=np.zeros_like(seg)
    a=np.unique(seg)
    pool=Pool(processes=nb_proc)
    count=0
    mapping=[]
    for i in a:
        tmp=seg[bboxes[i-1]]==i
        size_cell=np.sum(tmp)
        if size_cell>min_size_cell:
            count+=1
            mapping.append((tmp, \
                            max_size_cell, size_cell, \
                            iterations, out_iterations, \
                            bboxes[i-1], i))
    outputs=pool.map(create_seed, mapping)
    pool.close()
    pool.terminate()    
    for seed, num, bb in outputs:
        seeds[bb][seed]=num                         
    return SpatialImage(seeds, voxelsize=seg.voxelsize)



def cell_propagation(parameters):
    """
    Return the seeds in seeds_not_prop stricly included in cell c in seg_c
    seg_c : segmented image (SpatialImage)
    c : label of the cell to process
    seeds_not_prop : image of seeds (SpatialImage)
    """
    seg_c, c, seeds_not_prop=parameters
    if len(np.unique(seg_c))!=2: # DON'T MESS WITH THE SEEDS ! YOU NEED ONE AND ONLY ONE !!
        return
    seg_out=None
    labels=list(np.unique(seeds_not_prop[seg_c==c]))
    #if 0 in labels:
    labels.remove(0)#TODO Check if 0 is inside labels list 
    final_labels=[]
    for l in labels:
        if (seg_c[seeds_not_prop==l]==c).all():
            final_labels.append(l)
    nb=len(final_labels)
    return seg_out, nb, final_labels, c

def extract_seeds(seg_c, c, path_seeds_not_prop=None, bb=None, accept_3_seeds=False):
    """
    Return the seeds from seeds_not_prop stricly included in cell c from seg_c (the labels of the seeds go from 1 to 3)
    seg_c : segmented image (SpatialImage)
    c : label of the cell to process
    seeds_not_prop : image of seeds (can be the path to the image or the SpatialImage)
    bb : if seeds_not_prop is a path then bb is the bounding box of c in seeds_not_prop
    accept_3_seeds : True if 3 seeds can be accepted as a possible choice
    """
    if type(path_seeds_not_prop)!=SpatialImage:
        seeds_not_prop_out=imread(path_seeds_not_prop)
        seeds_not_prop=seeds_not_prop_out[bb]
    else: ## Then path_seeds_not_prop is the actual image we want to work with
        from copy import deepcopy
        seeds_not_prop=deepcopy(path_seeds_not_prop)
    labels=list(np.unique(seeds_not_prop[seg_c==c]))
    labels.remove(0)
    final_labels=[]
    for l in labels:
        if (seg_c[seeds_not_prop==l]==c).all():
            final_labels.append(l)    
    if len(final_labels)==1:
        return (1, (seeds_not_prop==final_labels[0]).astype(np.uint8))
    elif len(final_labels)==2:
        return (2, ((seeds_not_prop==final_labels[0]) + 
                    2*(seeds_not_prop==final_labels[1])).astype(np.uint8))
    elif len(final_labels)==3 and not accept_3_seeds:# "too much seeds in the second extraction"
        return (3, ((seeds_not_prop==final_labels[0]) + 
                    2*(seeds_not_prop==final_labels[1])).astype(np.uint8))
    elif len(final_labels)==3 and accept_3_seeds: #"accept 3 seeds !"
        return (3, ((seeds_not_prop==final_labels[0]) + 
                    2*(seeds_not_prop==final_labels[1]) +
                    3*(seeds_not_prop==final_labels[2])).astype(np.uint8))

def __slices_dilation(slices, maximum=[np.inf, np.inf, np.inf]):
    return tuple([slice(max(0, s.start-1), min(s.stop+1, maximum[i])) for i, s in enumerate(slices)])

def slices_dilation(slices, maximum=[np.inf, np.inf, np.inf], iterations=1):
    for i in range(iterations):
        slices=__slices_dilation(slices, maximum)
    return slices    


def to_u8(im, lt=0):
    """
     Return a SpatialImage in unsigned int
    im : SpatialImage
    lt : if the smallest value in the intensity image can be "predicted"
    """
    from copy import deepcopy
    imcp=deepcopy(im)
    tmp=imcp[:,:,imcp.shape[2]/3]
    fper=np.percentile(tmp[tmp>=lt], 1)
    nper=np.percentile(tmp[tmp>=lt], 99)
    imcp[imcp<fper]=fper
    imcp[imcp>nper]=nper
    #im-=fper
    np.subtract(imcp, fper, out=imcp, casting='unsafe')
    return SpatialImage(np.uint8(np.linspace(0, 255, nper-fper+1), casting='unsafe')[imcp], voxelsize=imcp.voxelsize)


def get_seeds(seg, h_min_min,h_min_max, sigma, cells, fused_file, path_h_min, bounding_boxes, nb_proc=26, verbose=False):
    """
    Return the number of seeds found for each cell in seg for different h_min values (from h_min_max down to 1)
    seg : Segmented image (SpatialImage)
    h_min_max : starting maximum value of h_min
    sigma : sigma of the gaussian smoothing (in voxels)
    cells : cells contained in seg
    fused_file : path (?) towards the fused image on which to perform the local minima detection
    path_h_min : format of h minima file names
    bounding_boxes : bounding boxes of the cells in seg (to fasten the computation)
    verbose : verbose mode (False or True)
    """
    from multiprocessing import Pool
    nb_cells={}
    treated=[]
    parameters={}
    mask=None
    temp_path_h_min=path_h_min.replace('$HMIN',str(h_min_max))
    if not os.path.exists(temp_path_h_min):
        seeds_not_prop, mask=obsolete_find_local_minima(temp_path_h_min, fused_file, h_min_max, sigma=sigma, verbose=verbose)
    else:
        seeds_not_prop=imread(temp_path_h_min)

    h_min=h_min_max
    #
    #
    #
    tmp_nb=[]
    checking=True
    while (checking):
        mapping=[]
        tmp_nb=[]
        for c in cells:
            if not c in treated:
                bb=slices_dilation(bounding_boxes[c], maximum=seg.shape, iterations=2)
                seg_c=np.ones_like(seg[bb])
                seg_c[seg[bb]==c]=c
                mapping.append((seg_c, c, seeds_not_prop[bb]))

        pool=Pool(processes=nb_proc)
        outputs=pool.map(cell_propagation, mapping)
        pool.close()
        pool.terminate()
        for seg_c_p, nb, labels, c in outputs:
            tmp_nb.append(nb)
            nb_cells.setdefault(c, []).append(nb)
            parameters.setdefault(c, []).append([h_min, sigma])

        h_min-=2
        checking=h_min>=h_min_min and (((np.array(tmp_nb)<=2) & (np.array(tmp_nb)!=0)).any() or tmp_nb==[])
        if checking :
            temp_path_h_min=path_h_min.replace('$HMIN',str(h_min))
            if not os.path.exists(temp_path_h_min):
                seeds_not_prop, mask=obsolete_find_local_minima(temp_path_h_min,fused_file, h_min, mask=mask, sigma=sigma, verbose=verbose)
            else:
                seeds_not_prop=imread(temp_path_h_min)
            if seeds_not_prop is None:
                checking=False
    return nb_cells, parameters


def get_back_parameters(nb_cells, parameters, lin_tree, cells,Thau=25):
    """
    Return the correct h-minima value for each cell
    nb_cells : { cell: [#seeds, ] }: dict, key: cell, values: list of #seeds
    parameters : { cell: [[h_min, sigma], ]}: dict matching nb_cells, key: cell, values: list of parameters
    """
    lin_tree_back={ v:k for k, val in lin_tree.iteritems() for v in val }
    right_parameters={}
    cells_with_no_seed=[]
    ## 2 plateau size vs noise ##
    for c, s in nb_cells.iteritems():
        nb_2=np.sum(np.array(s)==2)
        nb_3=np.sum(np.array(s)>=2)
        score=nb_2*nb_3
        if (s.count(1) or s.count(2))!=0:
            if score>=Thau:
                h, sigma=parameters[c][np.where(np.array(s)==2)[0][0]]
                nb_final=2
            elif s.count(1)!=0:
                h, sigma=parameters[c][np.where(np.array(s)==1)[0][0]]
                nb_final=1
            else:
                h, sigma=parameters[c][np.where(np.array(s)==2)[0][0]]
                nb_final=2
            right_parameters[c]=[h, sigma, nb_final]
        elif s.count(3)!=0:
            h, sigma=parameters[c][s.index(3)]
            right_parameters[c]=[h, sigma, 3]
        else:
            cells_with_no_seed.append(c)
            right_parameters[c]=[0, 0, 0]
    return right_parameters, cells_with_no_seed


def get_seeds_from_optimized_parameters(t, seg, cells, cells_with_no_seed, right_parameters,delta_t, bounding_boxes, im_ref, seeds, parameters, h_min_max, path_h_min, sigma,Volum_Min_No_Seed=100, verbose=False):
    """
    Return the seed image from the locally parametrized h-minima operator
    t : time
    seg : propagated segmentation (seg at t deformed on t+dt)
    cells : list of cells in seg
    cells_with_no_seed : list of cells with no correct parameters
    right_parameters : dict of the correct parameters for every cells
    delta_t : dt
    bounding_boxes : bounding boxes of the cells in seg (to fasten the computation)
    im_ref : Intensity image at time t+dt (on which to permorm the watershed)
    seeds : Propagated seeds from segmentation at time t (when no correct parameters were found)
    parameters : ?
    h_min_max : starting maximum value of h_min
    sigma : sigma of the gaussian smoothing (in voxels)
    path_h_min : format of h minima file names
    """
    seeds_from_opt_h=np.zeros_like(seg, dtype=np.uint16)
    label_max=2
    corres={}
    divided_cells=[]
    h_min_information={}
    sigma_information={}
    sigma_done=[]
    h_min_done=[]
    seeds_images={}
    for c in cells:
    	print 'get_seeds_from_optimized_parameters on '+str(c)
        if c in cells_with_no_seed:
            continue
        if not seeds_images.has_key((right_parameters[c][0], right_parameters[c][1])):
            path_seeds_not_prop=path_h_min.replace('$HMIN',str(right_parameters[c][0])).replace('$SIGMA',str(right_parameters[c][1]));
            seeds_images[(right_parameters[c][0], right_parameters[c][1])]=imread(path_seeds_not_prop)
        bb=slices_dilation(bounding_boxes[c], maximum=seg.shape, iterations=2)
        seg_c=np.ones_like(seg[bb])
        seg_c[seg[bb]==c]=c
        seeds_ex=seeds_images[(right_parameters[c][0], right_parameters[c][1])][bb]
        nb, seeds_c=extract_seeds(seg_c, c, seeds_ex)
        if nb==1:
            corres[c]=[label_max]
            h_min_information[(t+delta_t)*10**4+label_max]=right_parameters[c][0]
            sigma_information[(t+delta_t)*10**4+label_max]=right_parameters[c][1]
            seeds_from_opt_h[bb]+=seeds_c*label_max
            label_max+=1
        elif nb==2:
            corres[c]=[label_max, label_max+1]
            divided_cells.append((label_max, label_max+1))
            seeds_from_opt_h[bb][seeds_c==1]=label_max
            h_min_information[(t+delta_t)*10**4+label_max]=right_parameters[c][0]
            sigma_information[(t+delta_t)*10**4+label_max]=right_parameters[c][1]
            label_max+=1
            seeds_from_opt_h[bb][seeds_c==2]=label_max
            h_min_information[(t+delta_t)*10**4+label_max]=right_parameters[c][0]
            sigma_information[(t+delta_t)*10**4+label_max]=right_parameters[c][1]
            label_max+=1
        elif nb==3:
            corres[c]=[label_max, label_max+1]
            divided_cells.append((label_max, label_max+1))
            seeds_from_opt_h[bb][seeds_c==1]=label_max
            h_min_information[(t+delta_t)*10**4+label_max]=right_parameters[c][0]
            sigma_information[(t+delta_t)*10**4+label_max]=right_parameters[c][1]
            label_max+=1
            seeds_from_opt_h[bb][seeds_c==2]=label_max
            h_min_information[(t+delta_t)*10**4+label_max]=right_parameters[c][0]
            sigma_information[(t+delta_t)*10**4+label_max]=right_parameters[c][1]
            label_max+=1

	print 'Create Background seed'
    c=1
    seg_c=np.ones_like(seg)
    seg_c[seg!=c]=0
    sigma_out=sigma
    key_min = (h_min_max, sigma_out)
    for k in seeds_images.iterkeys():
    	if k[0]<key_min[0]:
    		key_min = k

    print 'Cell propagation'
    seeds_not_prop=seeds_images[key_min]
    parameters=(seg_c, c,seeds_not_prop)
    seg_c_p, nb, labels, c=cell_propagation(parameters)
    corres[1]=[]
    exterior_corres=[]
    for l in labels:
        seeds_from_opt_h=seeds_from_opt_h.astype(np.uint16)
        exterior_corres.append(label_max)
        seeds_from_opt_h[seeds_not_prop==l]=label_max
        label_max+=1

    print 'Cells with not Seed'
    for c in cells_with_no_seed:
        if np.sum(seg==c)>Volum_Min_No_Seed:
            seeds_from_opt_h[seeds==c]=label_max
            h_min_information[(t+delta_t)*10**4+label_max]=right_parameters[c][0]
            corres[c]=[label_max]
            label_max+=1

    print 'Watershed '
    seg_from_opt_h=obsolete_watershed(SpatialImage(seeds_from_opt_h, voxelsize=seeds_from_opt_h.voxelsize), im_ref, temporary_folder=os.path.dirname(path_h_min), verbose=verbose)
    for l in exterior_corres:
        seg_from_opt_h[seg_from_opt_h==l]=1
    corres[1]=[1]

    return seeds_from_opt_h, seg_from_opt_h, corres, exterior_corres, h_min_information, sigma_information, divided_cells, label_max    


def perform_ac(parameters):
    """
    Return the shape resulting of morphosnake operation on image I using image S as an initialisation
    m : label of the cell to work on
    daughters : list of the daughters of cell m (to keep track when working in parallel)
    bb : bounding boxe of m
    I : intensity image to perform active contours on (SpatialImage)
    S : segmented image to perform active contours from (SpatialImage, must contain the label m)
    """


    m, daughters, bb, I, S,MorphosnakeIterations,NIterations,DeltaVoxels=parameters
    import os
    from scipy import ndimage as nd
    import morphsnakes
    cell_num=m
    Sb=nd.binary_erosion(S!=cell_num, iterations=MorphosnakeIterations, border_value=1)#[:,:,sl]
    image_input='tmp_'+str(cell_num)+'.inr'
    gradient_output='tmp_out_'+str(cell_num)+'.inr'
    imsave(image_input, I)
    obsolete_gradient_norm(image_input,gradient_output)
    gI = imread(gradient_output)
    os.system('rm -f '+image_input+' '+gradient_output)
    gI=1./np.sqrt(1+100*gI)
	
	
    macwe = morphsnakes.MorphGAC(gI, smoothing=3, threshold=1, balloon=1)
    macwe.levelset = Sb
    bef=np.ones_like(Sb)
    from copy import deepcopy
    for i in xrange(NIterations):
        beff=deepcopy(bef)
        bef=deepcopy(macwe.levelset)
        macwe.step()
        if np.sum(bef!=macwe.levelset)<DeltaVoxels or np.sum(beff!=macwe.levelset)<DeltaVoxels:
            break
    out=macwe.levelset
    tmp=nd.binary_fill_holes(out)
    cell_out=(out.astype(np.bool) ^ tmp)
    return m, daughters, bb, cell_out


def volume_checking(t,delta_t,seg, seeds_from_opt_h, seg_from_opt_h, corres, divided_cells, bounding_boxes, right_parameters, im_ref, im_ref16, seeds, nb_cells, 
    label_max, exterior_corres, parameters, h_min_information, sigma_information, segmentation_file_ref, segmentation_file_trsf, path_h_min, volumes_t_1, 
    nb_proc=26,Thau= 25,MinVolume=1000,VolumeRatioBigger=0.5,VolumeRatioSmaller=0.1,MorphosnakeIterations=10,NIterations=200 ,DeltaVoxels=10**3):
    """
    Return corrected final segmentation based on conservation of volume in time
    seg : propagated segmentation (seg at t deformed on t+dt) (SpatialImage)
    seeds_from_opt_h : optimized seeds (SpatialImage)
    seg_from_opt_h : segmented image from seeds_from_opt_h (SpatialImage)
    corres : mapping of cells at time t to cells at t+dt in seg_from_opt_h
    divided_cells : list of cells that have divided between t and t+dt
    bounding_boxes : bounding boxes of the cells in seg (to fasten the computation)
    right_parameters : list of parameters used to create seeds_from_opt_h
    im_ref : image to segment at time t+dt 8 bits (SpatialImage)
    im_ref16 : image to segment at time t+dt in 16 bits (SpatialImage)
    seeds : Propagated seeds from segmentation at time t
    nb_cells : { cell: [#seeds, ] }: dict, key: cell, values: list of #seeds
    label_max : maximum label in seg_from_opt_h
    exterior_corres : list of cells that have been corrected for issue in exterior
    parameters : { cell: [[h_min, sigma], ]}: dict matching nb_cells, key: cell, values: list of parameters
    h_min_information : { cell: h_min}: dict associating to each cells the h_min that allowed its segmentation
    sigma_information : { cell: sigma}: dict associating to each cells the sigma that allowed its segmentation
    segmentation_file_ref : path to the segmentation at time t
    segmentation_file_trsf : path to the segmentation at time t resampled at t+1 
    vf_file : path to the vector field that register t into t+dt
    path_h_min : format of h-minima files
    volumes_t_1 : cell volumes at t
    """

    # seg_origin : original segmentation (SpatialImage)

    seg_origin=imread(segmentation_file_ref)

    volumes_from_opt_h=compute_volumes(seg_from_opt_h)
    if volumes_t_1=={}:
        volumes=compute_volumes(seg_origin)
    else:
        volumes=volumes_t_1

    bigger=[]
    lower=[]
    to_look_at=[]
    too_little=[]
    for mother_c, sisters_c in corres.iteritems():
        if mother_c!=1:
            volume_ratio=1-(volumes[mother_c]/np.sum([volumes_from_opt_h.get(s, 1) for s in sisters_c]))
            if not (-VolumeRatioSmaller<volume_ratio<VolumeRatioSmaller):
                if (volume_ratio>0) and (volumes_from_opt_h.get(s, 1)!=1):
                    bigger.append((mother_c, sisters_c))
                elif volumes_from_opt_h.get(s, 1)!=1 :
                    lower.append((mother_c, sisters_c))
                if volume_ratio<-VolumeRatioBigger:
                    to_look_at.append(mother_c)
            else :
                for s in sisters_c:
                    if volumes_from_opt_h[s]<MinVolume:
                        too_little.append((mother_c, s))
                
    to_fuse_3=[]
    change_happen=False
    for c in to_look_at:
        s=nb_cells[c]
        nb_2=np.sum(np.array(s)==2)
        nb_3=np.sum(np.array(s)>=2)
        score=nb_2*nb_3
        if (s.count(1) or s.count(2))!=0:
            if score>=Thau:
                h, sigma=parameters[c][np.where(np.array(s)==2)[0][-1]]
                nb_final=2
            elif s.count(1)!=0:
                h, sigma=parameters[c][np.where(np.array(s)==1)[0][-1]]
                nb_final=1
            else:
                h, sigma=parameters[c][np.where(np.array(s)==2)[0][-1]]
                nb_final=2
            right_parameters[c]=[h, sigma, nb_final]

            if nb_final==1 and s.count(2)!=0:
                h, sigma=parameters[c][s.index(2)]
                path_seeds_not_prop=path_h_min.replace('$HMIN',str(h)).replace('$SIGMA',str(sigma));
                bb=slices_dilation(bounding_boxes[c], maximum=seg.shape, iterations=2)
                seg_c=np.ones_like(seg[bb])
                seg_c[seg[bb]==c]=c
                nb, seeds_c=extract_seeds(seg_c, c, path_seeds_not_prop, bb)

                if nb==2 and (seg_from_opt_h[bb][seeds_c!=0]==0).any(): #If we can found 2 seeds and one is not in new calculated cell
                    change_happen=True
                    seeds_from_opt_h[seeds_from_opt_h==corres[c][0]]=0
                    corres[c]=[label_max, label_max+1]
                    divided_cells.append((label_max, label_max+1))
                    seeds_from_opt_h[bb][seeds_c==1]=label_max
                    h_min_information[(t+delta_t)*10**4+label_max]=h
                    sigma_information[(t+delta_t)*10**4+label_max]=sigma
                    label_max+=1
                    seeds_from_opt_h[bb][seeds_c==2]=label_max
                    h_min_information[(t+delta_t)*10**4+label_max]=h
                    sigma_information[(t+delta_t)*10**4+label_max]=sigma
                    label_max+=1
            if (nb_final==1 or nb_final==2) and (np.array(s)>2).any():
                h, sigma=parameters[c][-1]
                path_seeds_not_prop=path_h_min.replace('$HMIN',str(h)).replace('$SIGMA',str(sigma));
                seeds_image=imread(path_seeds_not_prop)
                bb=slices_dilation(bounding_boxes[c], maximum=seg.shape, iterations=2)
                seg_c=np.zeros_like(seg_from_opt_h[bb])
                for daughter in corres[c]:
                    seg_c[seg_from_opt_h[bb]==daughter]=1
                seeds_c=np.zeros_like(seg_from_opt_h[bb])
                seeds_c[(seg_c==1) & (seeds_image[bb]!=0)]=1
                seeds_c[(seg[bb]==c) & (seg_c!=1) & (seeds_image[bb]!=0)]=2
                if 2 in seeds_c:
                    change_happen=True
                    for daughter in corres[c]:
                        seeds_from_opt_h[seeds_from_opt_h==daughter]=0
                    corres[c]=[label_max, label_max+1]
                    divided_cells.append((label_max, label_max+1))
                    seeds_from_opt_h[bb][seeds_c==1]=label_max
                    h_min_information[(t+delta_t)*10**4+label_max]=h
                    sigma_information[(t+delta_t)*10**4+label_max]=sigma
                    label_max+=1
                    seeds_from_opt_h[bb][seeds_c==2]=label_max
                    h_min_information[(t+delta_t)*10**4+label_max]=h
                    sigma_information[(t+delta_t)*10**4+label_max]=sigma
                    label_max+=1
            elif nb_final==1:
                change_happen=True
                seeds_from_opt_h[seeds_from_opt_h==corres[c][0]]=0
                seeds_from_opt_h[seeds==c]=corres[c][0]
                label_max+=1
        elif s.count(3)!=0:
            h, sigma=parameters[c][s.index(3)]
            path_seeds_not_prop=path_h_min.replace('$HMIN',str(h)).replace('$SIGMA',str(sigma));
            bb=slices_dilation(bounding_boxes[c], maximum=seg.shape, iterations=2)
            seg_c=np.ones_like(seg[bb])
            seg_c[seg[bb]==c]=c
            nb, seeds_c=extract_seeds(seg_c, c, path_seeds_not_prop, bb, accept_3_seeds=True)
            change_happen=True
            #addition to correct 0-boolean error when len(corres[c])>1
            for ci in range(len(corres[c])):
                seeds_from_opt_h[seeds_from_opt_h==corres[c][ci]]=0
            #seeds_from_opt_h[seeds_from_opt_h==corres[c]]=0
            divided_cells.append((label_max, label_max+1))
            seeds_from_opt_h[bb][seeds_c==1]=label_max
            h_min_information[(t+delta_t)*10**4+label_max]=h
            sigma_information[(t+delta_t)*10**4+label_max]=sigma
            label_max+=1
            seeds_from_opt_h[bb][seeds_c==2]=label_max
            h_min_information[(t+delta_t)*10**4+label_max]=h
            sigma_information[(t+delta_t)*10**4+label_max]=sigma
            label_max+=1
            seeds_from_opt_h[bb][seeds_c==3]=label_max
            h_min_information[(t+delta_t)*10**4+label_max]=h
            sigma_information[(t+delta_t)*10**4+label_max]=sigma
            label_max+=1
            to_fuse_3.append([c, (label_max-1, label_max-2, label_max-3)])



    if too_little!=[]:
        for c in too_little:
            #for d in corres[c]:
            seeds_from_opt_h[seeds_from_opt_h==c[1]]=0
            tmp=corres[c[0]]
            tmp.remove(c[1])
            if tmp==[]:
                corres.pop(c[0])
            else:
                corres[c[0]]=tmp
            change_happen=True

    if change_happen:
        seg_from_opt_h=obsolete_watershed(SpatialImage(seeds_from_opt_h,voxelsize=seeds_from_opt_h.voxelsize), im_ref, temporary_folder=os.path.dirname(path_h_min))
        for l in exterior_corres:
            seg_from_opt_h[seg_from_opt_h==l]=1  
            
        volumes_from_opt_h=compute_volumes(seg_from_opt_h)

    lower=[]
    for mother_c, sisters_c in corres.iteritems():
        if mother_c!=1:
            volume_ratio=1-(volumes[mother_c]/np.sum([volumes_from_opt_h.get(s, 1) for s in sisters_c]))
            if not (-.1<volume_ratio<.1):
                if (volume_ratio<0) and volumes_from_opt_h.get(s, 1)!=1 :
                    lower.append((mother_c, sisters_c))
    
    exterior_correction=[]
    if lower!=[]:
        from copy import deepcopy
        #tmp=apply_trsf(segmentation_file_ref, vf_file, nearest=True, lazy=False)
        tmp=imread(segmentation_file_trsf)
        old_bb=nd.find_objects(tmp)
        for mother_c, sisters_c in lower:
            cell_before=tmp[old_bb[mother_c-1]]==mother_c
            cell_after=np.zeros_like(cell_before)
            for c in sisters_c:
                cell_after+=seg_from_opt_h[old_bb[mother_c-1]]==c
            lost=seg_from_opt_h[old_bb[mother_c-1]][cell_after^cell_before]
            max_share=0
            share_lab=0
            size={}
            for v in np.unique(lost):
                size[v]=np.sum(lost==v)
                if np.sum(lost==v)>max_share:
                    max_share=np.sum(lost==v)
                    share_lab=v
            if share_lab==1 and 1 in tmp[old_bb[mother_c-1]]:
                exterior_correction.append((mother_c, sisters_c))
        from multiprocessing import Pool
        pool=Pool(processes=nb_proc)
        mapping=[]
        for m, daughters in exterior_correction:
            bb=slices_dilation(old_bb[m-1], maximum=im_ref.shape, iterations=15)
            im_ref_tmp=deepcopy(im_ref16[bb])
            seg_ref_tmp=deepcopy(tmp[bb])
            mapping.append((m, daughters, bb, im_ref_tmp, seg_ref_tmp,MorphosnakeIterations,NIterations,DeltaVoxels))
        outputs=pool.map(perform_ac, mapping)
        pool.close()
        pool.terminate()
        for m, daughters, bb, cell_out in outputs:
            seg_from_opt_h[bb][seg_from_opt_h[bb]==1 & cell_out]=daughters[0]
            if len(daughters)==2:
                seg_from_opt_h[bb][seg_from_opt_h[bb]==daughters[1]]=daughters[0]
                if tuple(daughters) in divided_cells:
                    divided_cells.remove(tuple(daughters))
            corres[m]=[daughters[0]]
    for c, tf in to_fuse_3:
        bb=slices_dilation(bounding_boxes[c], maximum=seg.shape, iterations=2)
        seg_c=np.ones_like(seg_from_opt_h[bb])
        seg_c[seg_from_opt_h[bb]==tf[0]]=tf[0]
        seg_c[seg_from_opt_h[bb]==tf[1]]=tf[1]
        seg_c[seg_from_opt_h[bb]==tf[2]]=tf[2]
        v1=np.sum(seg_c==tf[0])
        v2=np.sum(seg_c==tf[1])
        v3=np.sum(seg_c==tf[2])
        vol_cells_to_f=[v1, v2, v3]
        cell_to_f=np.argmin(vol_cells_to_f)
        tmp=nd.binary_dilation(seg_c==tf[cell_to_f])
        p1=tf[np.argsort(vol_cells_to_f)[1]]
        p2=tf[np.argsort(vol_cells_to_f)[2]]
        im_tmp=np.zeros_like(seg_c)
        im_tmp[seg_c==p1]=p1
        im_tmp[seg_c==p2]=p2
        im_tmp[tmp==False]=0
        p1_share=np.sum(im_tmp==p1)
        p2_share=np.sum(im_tmp==p2)
        if p1_share>p2_share:
            seg_from_opt_h[seg_from_opt_h==tf[cell_to_f]]=p1
        else:
            seg_from_opt_h[seg_from_opt_h==tf[cell_to_f]]=p2
        corres[c]=[p1, p2]
        divided_cells.append((p1, p2))

    return seg_from_opt_h, bigger, lower, to_look_at, too_little, corres, exterior_correction


def outer_correction(seg_from_opt_h, exterior_correction,segmentation_file_ref,RadiusOpening=20):
    """
    Return an eroded segmentation correcting for potential errors in the morphsnake
    seg_from_opt_h : segmentated image (SpatialImage)
    exterior_correction : list of cells that have been corrected using morphsnake algorithm
    """
    if exterior_correction!=[]:
    	image_input=segmentation_file_ref.replace('.inr','.seg_from_opt_h.inr')
        imsave(image_input, SpatialImage(seg_from_opt_h!=1, voxelsize=seg_from_opt_h.voxelsize).astype(np.uint8))
        image_output=segmentation_file_ref.replace('.inr','.seg_out_h.inr')
        obsolete_morpho(image_input,image_output,' -ope -R '+str(RadiusOpening))
        opened=imread(image_output)
        cells_to_correct=[i for j in exterior_correction for i in j[1]]
        os.system('rm -f '+image_input+' '+image_output)
        to_remove=opened^(seg_from_opt_h>1)
        
        for c in cells_to_correct:
            seg_from_opt_h[((seg_from_opt_h==c) & to_remove).astype(np.bool)]=1
    return seg_from_opt_h


def segmentation_propagation_seeds_init_and_deform(t, segmentation_ref, fused_file, seeds_file, vf_file, delta_t, verbose=False):
    """
    Steps 2 to 3 of segmentation propagation:
    create seeds from reference segmentation, then resample it by transformation application
    -> generation of seeds_file containing the image of seeds which will be used for segmentation propagation at t+delta_t
    """
    print 'Create The Seeds from '+str(t)

    seeds_ref=create_seeds(segmentation_ref, max_size_cell=np.inf)
    imsave(seeds_file, SpatialImage(seeds_ref, voxelsize=seeds_ref.voxelsize))
    
    print 'Deform Seeds with vector fields from '+str(t)+' to '+str(t+delta_t)
    apply_trsf(seeds_file, vf_file , path_output=seeds_file, template=fused_file,nearest=True, lazy=True, verbose=verbose)



def segmentation_propagation_from_seeds(t, segmentation_file_ref, fused_file,  fused_file_u8 , seeds_file,path_seg_trsf, path_h_min, h_min_min,h_min_max, sigma, lin_tree_information, delta_t, nb_proc,
    RadiusOpening=20,Thau=25,MinVolume=1000,VolumeRatioBigger=0.5,VolumeRatioSmaller=0.1,MorphosnakeIterations=10,NIterations=200,DeltaVoxels=10**3,Volum_Min_No_Seed=100, delSeedsASAP=True, 
    verbose=False):
    """
    Steps 4 to 9 of segmentation propagation:
    - initial watershed
    - computation of h-minima (get_seeds method)
    - optimal h selection for each cell (get_back_parameters method)
    - build a seeds image from previous information and new segmentation by watershed (get_seeds_from_optimized_parameters)
    - morphosnake if needed (called from volume_checking method)
    - last corrections with a morphological opening (outer_correction method)
    Returns seg_from_opt_h, lin_tree_information
        seg_from_opt_h : SpatialImage of the segmentation at t+delta_t
        lin_tree_information : updated lineage tree 
    """
    from copy import deepcopy
    lin_tree=lin_tree_information.get('lin_tree', {})
    tmp=lin_tree_information.get('volumes_information', {})
    volumes_t_1={k%10**4: v for k, v in tmp.iteritems() if k/10**4 == t}
    h_min_information={}
    

    print 'Perform watershed with the seeds from method "segmentation_propagation_seeds_init_and_deform"'
    im_fused=imread(fused_file)
    im_fused_8=imread(fused_file_u8)


    #
    # segmentation par watershed a l'aide des graines extraites de la segmentation precedente
    #
    segmentation=obsolete_watershed(seeds_file, im_fused_8, temporary_folder=os.path.dirname(path_seg_trsf), verbose=verbose)
    seeds=imread(seeds_file)
    if delSeedsASAP:
        cmd='rm %s'%seeds_file
        if verbose:
            print cmd
        os.system(cmd)

    #
    # bounding boxes
    #
    cells=list(np.unique(segmentation))
    cells.remove(1)
    bounding_boxes=dict(zip(range(1, max(cells)+1), nd.find_objects(segmentation)))
    treated=[]

    print 'Estimation of the local h-minimas at '+str(t+delta_t)
    nb_cells, parameters=get_seeds(segmentation, h_min_min,h_min_max, sigma, cells, fused_file,
                                   path_h_min.replace('$SIGMA',str(sigma)), bounding_boxes, nb_proc=nb_proc, verbose=verbose)

    #
    #
    #
    right_parameters, cells_with_no_seed=get_back_parameters(nb_cells, parameters, lin_tree, cells,Thau=Thau)
    
    print 'Applying volume correction '+str(t+delta_t)
    seeds_from_opt_h, seg_from_opt_h, corres, exterior_corres, h_min_information, sigma_information, divided_cells, label_max = get_seeds_from_optimized_parameters(t, segmentation, cells, cells_with_no_seed, 
        right_parameters, delta_t, bounding_boxes, im_fused_8, seeds, parameters, h_min_max, path_h_min, sigma,Volum_Min_No_Seed=Volum_Min_No_Seed, verbose=verbose)
    
    print 'Perform volume checking '+str(t+delta_t)
    seg_from_opt_h, bigger, lower, to_look_at, too_little, corres, exterior_correction = volume_checking(t,delta_t,segmentation, seeds_from_opt_h, seg_from_opt_h, corres, divided_cells, bounding_boxes, right_parameters, 
        im_fused_8, im_fused, seeds, nb_cells, label_max, exterior_corres, parameters, h_min_information, sigma_information, segmentation_file_ref, path_seg_trsf, path_h_min, volumes_t_1, 
        nb_proc=nb_proc,Thau=Thau, MinVolume=MinVolume,VolumeRatioBigger=VolumeRatioBigger,VolumeRatioSmaller=VolumeRatioSmaller,MorphosnakeIterations=MorphosnakeIterations,NIterations=NIterations ,DeltaVoxels=DeltaVoxels)

    print 'Perform Outer Correction '+str(t+delta_t)
    seg_from_opt_h = outer_correction(seg_from_opt_h, exterior_correction,segmentation_file_ref,RadiusOpening=RadiusOpening)

    print 'Compute Volumes'+str(t+delta_t)
    volumes=compute_volumes(seg_from_opt_h)
    volumes_information={}
    for k, v in volumes.iteritems():
        volumes_information[(t+delta_t)*10**4+k]=v
    for m, d in corres.iteritems():
        if m!=1:
            daughters=[]
            for c in d:
                if c in volumes:
                    daughters.append(c+(t+delta_t)*10**4)
                else:
                    print str(c) +' is not segmented'
            if len(daughters)>0:
                lin_tree[m+t*10**4]=daughters
    lin_tree_information['lin_tree']=lin_tree
    lin_tree_information.setdefault('volumes_information', {}).update(volumes_information)
    lin_tree_information.setdefault('h_mins_information', {}).update(h_min_information)
    lin_tree_information.setdefault('sigmas_information', {}).update(sigma_information)

    return seg_from_opt_h, lin_tree_information


def segmentation_propagation(t, fused_file_ref, segmentation_file_ref, fused_file , seeds_file,vf_file, path_h_min,
                             h_min_min, h_min_max, sigma, lin_tree_information, delta_t, nb_proc,
    membrane_reconstruction_method=None, fusion_u8_method=0, flag_hybridation=False, 
    RadiusOpening=20,Thau=25,MinVolume=1000,VolumeRatioBigger=0.5,VolumeRatioSmaller=0.1,MorphosnakeIterations=10,NIterations=200,DeltaVoxels=10**3,Volum_Min_No_Seed=100, 
    rayon_dil=3.6, sigma_membrane=0.9, manual=False, manual_sigma=7, hard_thresholding=False, hard_threshold=1.0, sensitivity=0.99, sigma_TV=3.6, sigma_LF=0.9, sample=0.2, 
    keep_membrane=False, keep_all=False,  path_u8_images=None, nb_proc_ACE=7, 
    min_percentile=0.01, max_percentile=0.99, min_method='cellinterior', max_method='cellborder', sigma_hybridation=5.0, 
    verbose=False, keepTemporaryFiles=False ):
    '''
    Return the propagated segmentation at time t+dt and the updated lineage tree and cell informations
    t : time t
    fused_file_ref : path format to fused images
    segmentation_file_ref : path format to segmentated seeds_images
    fused_file : fused image at t+dt
    vf_file : path format to transformation
    path_h_min : path format to h-minima files
    h_min_max : maximum value of the h-min value for h-minima operator
    sigma : sigma value in voxels for gaussian filtering
    lin_tree_information : dictionary containing the lineage tree dictionary, volume information, h_min information and sigma information for every cells
    delta_t : value of dt (in number of time point)
    nb_proc : number maximum of processors to allocate

    # Modules choice

    membrane_reconstruction_method : if not set or set to 0, the input fused_file is not processed for membrane structures enhancement.
                                     if set to 1, the GLACE reconstruction method is going to be called
                                     if set to 2, the GACE reconstruction method is going to be called

    fusion_u8_method : select method to convert fused_file into a 8 bits images for the segmentation propagation. 
                       if set to 0 (default), calling the historical "to_u8" method
                       if set to 1, calling the mc_adhocFuse function which enhances the fused image while converting it to u8 
                       knowing the semgnetation propagation from previous time point

    flag_hybridation : if set to True and if the membrane_reconstruction_method parameter is provided and not equal to 0, 
                       then the reconstructed gray level image
                       used for semgentation_propragation_from_seeds is goind to be ahybridation between the original image
                       fused_file and the result of image reconstruction by the specified method.
    
    path_u8_images : default is None. If provided, saves a copy of the u8 image used for watershed process at the specified path.



    # Glace Parameters (if membrane_reconstruction_method is set to 1 or 2):
    # membrane_renforcement
    sigma_membrane=0.9  # membrane enhancement parameter (in real units, a
                        # priori 0.9 um is a good choice for data like 
                        # Patrick/Ralph/Aquila)
    # anisotropicHist /!\ critical step
    sensitivity=0.99    # membrane binarization parameter, /!\ if failure,
                        # one should enter in "manual" mode of the function
                        # anisotropicHist via activation of 'manual' option

    manual=False        # By default, this parameter is set to False. If 
                        # failure, (meaning that thresholds are very bad, 
                        # meaning that the binarized image is very bad),
                        # set this parameter to True and relaunch the 
                        # computation on the test image. If the method fails
                        # again, "play" with the value of manual_sigma... 
                        # and good luck.
    manual_sigma=15     # Axial histograms fitting initialization parameter 
                        # for the computation of membrane image binarization
                        # axial thresholds (this parameter is used iif 
                        # manual = True).
                        # One may need to test different values of 
                        # manual_sigma. We suggest to test values between 5 and
                        # 25 in case of initial failure. Good luck.

    hard_thresholding=False   # If the previous membrane threshold method 
                              # failed, one can force the thresholding with a
                              # "hard" threshold applied on the whole image. 
                              # To do so, this option must be set to True.
    hard_threshold=1.0        # If hard_thresholding = True, the enhanced 
                              # membranes image is thresholded using this 
                              # parameter (value 1 seems to be ok for 
                              # time-point t001 of Aquila embryo for example).

    # Tensor voting framework
    sigma_TV=3.6    # parameter which defines the voting scale for membrane
                    # structures propagation by tensor voting method (real
                    # coordinates). 
                    # This parameter shoud be set between 3 um (little cells)
                    # and 4.5 um(big gaps in the binarized membrane image)
    sigma_LF=0.9    # Smoothing parameter for reconstructed image (in real
                    # coordinates). It seems that the default value = 0.9 um
                    # is ok for classic use.
    sample=0.2      # Parameter for tensor voting computation speed 
                    # optimisation (do not touch if not bewared)
    rayon_dil=3.6   # dilatation ray for propagated ROI from time t to t+1
                    # (default: 3.6, in real coordinates) 

    nb_proc_ACE=7   # number of processors for ACE (7 is recommanded)

    '''
    segmentation_ref=imread(segmentation_file_ref);


    print 'Compute Vector Fields from '+str(t)+' to '+str(t+delta_t)
    obsolete_non_linear_registration(fused_file_ref,\
                        fused_file, \
                        vf_file.replace('.inr','_affine.inr'), \
                        vf_file.replace('.inr','_affine.trsf'),\
                        vf_file.replace('.inr','_vector.inr'),\
                        vf_file);

    cmd='rm -f '+vf_file.replace('.inr','_affine.inr')+' '+vf_file.replace('.inr','_affine.trsf')+' '+vf_file.replace('.inr','_vector.inrf')
    if verbose:
        print cmd
    os.system(cmd)

    segmentation_propagation_seeds_init_and_deform(t, segmentation_ref, fused_file, seeds_file, vf_file, delta_t, verbose=verbose)


    # graylevel image construction for segmentation propagation

    # defining temporary file paths 
    graylevel_file=vf_file.replace('.inr','_graylevel.inr')         # The first input gray level image for segmentation_propagation_from_seeds
    graylevel_file_u8=vf_file.replace('.inr','_graylevel_u8.inr')   # The second input gray level image for segmentation_propagation_from_seeds (must be a 8 bits image)
    fused_file_u8=vf_file.replace('.inr','_fuse_u8.inr')            # Temporary file
    path_seg_trsf=vf_file.replace('.inr','_seg_trsf.inr')           # Temporary file

    # segmentation propagation 
    apply_trsf(segmentation_file_ref, path_trsf=vf_file, path_output=path_seg_trsf, template=fused_file, nearest=True, verbose=verbose)

    # transformation file deletion
    if keepTemporaryFiles == False:
        cmd='rm -f '+vf_file
        if verbose:
            print cmd
        os.system(cmd)


    # fused file u8 vconversion if needed
    if flag_hybridation or not membrane_reconstruction_method:
        if  fusion_u8_method==1:
            obsolete_mc_adhocFuse(fused_file, path_seg_trsf, fused_file_u8, min_percentile=min_percentile, max_percentile=max_percentile,
                         min_method=min_method, max_method=max_method, sigma=sigma_hybridation, verbose=verbose)
        else:
            imsave(fused_file_u8, to_u8(imread(fused_file)))

    # Switch membrane_reconstruction_method
    if not membrane_reconstruction_method:
        obsolete_copy(fused_file_u8, graylevel_file_u8, verbose=verbose)
        obsolete_copy(fused_file, graylevel_file, verbose=verbose)
    if membrane_reconstruction_method == 1:
        # GLACE reconstruction 
        GLACE_from_resampled_segmentation(fused_file, path_seg_trsf, labels_of_interest='all', background=[0,1], 
        path_output=graylevel_file, rayon_dil=rayon_dil, 
        sigma_membrane=sigma_membrane, manual=manual, manual_sigma=manual_sigma, hard_thresholding=hard_thresholding, 
        hard_threshold=hard_threshold, sensitivity=sensitivity, sigma_TV=sigma_TV, sigma_LF=sigma_LF, sample=sample, 
        keep_membrane=keep_membrane, keep_all=keep_all,  nb_proc=nb_proc_ACE, verbose=verbose)   
    if membrane_reconstruction_method == 2:
        # GACE reconstruction
        out=GACE(fused_file, binary_input=False, path_output=graylevel_file, 
        sigma_membrane=sigma_membrane, manual=manual, manual_sigma=manual_sigma, hard_thresholding=hard_thresholding, 
        hard_threshold=hard_threshold, sensitivity=sensitivity, sigma_TV=sigma_TV, sigma_LF=sigma_LF, sample=sample, 
        keep_membrane=keep_membrane, keep_all=keep_all, verbose=verbose)



    # reconstructed image and fused image hybridation if needed
    if membrane_reconstruction_method:
        if flag_hybridation:
            obsolete_Arit(fused_file_u8, graylevel_file, graylevel_file, Mode='max', Type='-o 1', verbose=verbose)
            obsolete_copy(graylevel_file, graylevel_file_u8, verbose=verbose)

    # temporary images deletion
    if keepTemporaryFiles == False:
        if os.path.exists(fused_file_u8):
            cmd='rm -f '+fused_file_u8
            if verbose:
                print cmd
            os.system(cmd)

    # u8 image copy if asked
    if path_u8_images:
        obsolete_copy(graylevel_file_u8, path_u8_images, verbose=verbose)

    # segmentation propagation stuff from seeds
    seg_from_opt_h, lin_tree_information = segmentation_propagation_from_seeds(t, segmentation_file_ref, graylevel_file, graylevel_file_u8, seeds_file,path_seg_trsf, 
                                                                               path_h_min, h_min_min,h_min_max, sigma, lin_tree_information, delta_t, nb_proc,
                                                                               RadiusOpening=RadiusOpening,Thau=Thau,MinVolume=MinVolume,VolumeRatioBigger=VolumeRatioBigger,
                                                                               VolumeRatioSmaller=VolumeRatioSmaller,MorphosnakeIterations=MorphosnakeIterations,
                                                                               NIterations=NIterations,DeltaVoxels=DeltaVoxels,Volum_Min_No_Seed=Volum_Min_No_Seed, 
                                                                               delSeedsASAP=True, verbose=verbose)

    # temporary images deletion
    if keepTemporaryFiles == False:
        if os.path.exists(path_seg_trsf):
            cmd='rm -f '+path_seg_trsf
            if verbose:
                print cmd
            os.system(cmd)

        if os.path.exists(graylevel_file):
            cmd='rm -f '+graylevel_file
            if verbose:
                print cmd
            os.system(cmd)

        if os.path.exists(graylevel_file_u8):
            cmd='rm -f '+graylevel_file_u8
            if verbose:
                print cmd
            os.system(cmd)

    return seg_from_opt_h, lin_tree_information
