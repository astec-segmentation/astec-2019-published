import os, sys

assert os.path.isdir(os.path.join(os.path.dirname(__file__),"CommunFunctions"))
sys.path.append(os.path.join(os.path.dirname(__file__),"CommunFunctions"))

from scipy import ndimage as nd
import numpy as np
from copy import copy,deepcopy
from ImageHandling import imread, imsave, SpatialImage
import pickle as pkl
from lineage import read_lineage_tree,write_lineage_tree,timesNamed

def timeNamed(filename, time):
    time_point = ('00' + str(time))[-3:]  # Format time on 3 digit
    return filename.replace('$TIME', time_point)+".inr"


def fuse_cells(couple, lin_tree, bb, im, t):
    """
    Return an segmented image where the cells in couple have been fused
    couple : tuple of cells to fuse, couple[0] is fused to couple[1]
    lin_tree : lineage tree
    bb : bounding box of the couple of cells in ``couple``
    im : segmented image to fuse
    t : time point of im
    """
    inv_lin_tree={ v : k for k, values in lin_tree.iteritems() for v in values }
    give=couple[0]%10**4
    get=np.array(couple[1])%10**4
    possibles=np.array(lin_tree[inv_lin_tree[t*10**4+give]])
    possibles=possibles[possibles!=10**4*t+give]%10**4
    possibles=np.array([p for p in possibles if p in get])
    get=possibles
    if len(get)==1:
        final_get=get[0]
    else:
        if bb is None:
            bb=nd.find_objects(im)
        im_tmp=im[bb[give-1]]
        im_tmp_d=nd.binary_dilation(im_tmp==give)
        borders_b=im_tmp_d & (im_tmp!=give)
        borders=im_tmp[borders_b]
        count=[]
        for g in get:
            count.append([g, np.sum(borders==g)])
        count=np.array(count)
        final_get=count[np.argmax(count[:,1]), 0]
    return bb, (give, final_get)
    


def apply_cell_fusion(lin_tree,volumes, to_fuse, segmented_files,segmented_fused_files,begin, end, step):
    """
    Fuse all the cells from the information given by the lineage correction for a time series
    begin : starting point
    end : ending point
    step : dt
    path : path to the working folder
    folder : folder containing the segmentations (in path)
    path_lin_tree : path to the lineage tree
    to_fuse_file : path the the file containing the informations on the cells to fuse
    path_file_format : format of the segmented files
    """
    print 'Process cell fusion'

    inv_lin_tree={ v : k for k, values in lin_tree.iteritems() for v in values }

    for t in range(begin, end+1, step):
        current_state=deepcopy(lin_tree)
        print ' ->Cell fusion at '+str(t)
        bb=None
        time=('00'+str(t))[-3:]
        time_prev=('00'+str(t-1))[-3:]
        tf=to_fuse.get(t, '')
        ext_cor=[]
        segmented_file=timeNamed(segmented_files,t)
        im=imread(segmented_file)
        treated=[]
        mapping=[]
        if tf!='':
            if im is None:
                im=imread(file_name)
            if mapping==[]:
                mapping=np.array(range(np.max(im)+1), dtype=np.uint16)
            getting=[]
            for couple in tf:
                bb, trsf=fuse_cells(couple, lin_tree, bb, im, t)
                getting.append(trsf)
            getting={i:j for i,j in getting}
            roots=set(getting.keys()).difference(set(getting.values()))
            final_pairing={}
            for r in roots:
                n=r
                to_change=[]
                while getting.has_key(n):
                    to_change.append(n)
                    n=getting[n]
                    final=n
                for n in to_change:
                    getting[n]=final
            for give, get in getting.iteritems():
                inv_lin_tree={ v : k for k, values in lin_tree.iteritems() for v in values }
                if give<mapping.shape[0]:
                    mapping[give]=get
                tmp=lin_tree.get(10**4*t+get, [])
                tmp.extend(lin_tree.get(10**4*t+give, []))
                lin_tree[inv_lin_tree[10**4*t+get]].remove(10**4*t+give)
                lin_tree[10**4*t+get]=tmp
                lin_tree.pop(10**4*t+give, None)
                volumes[10**4*t+get]=volumes[10**4*t+get]+volumes[10**4*t+give]
        if mapping!=[]:
            im=SpatialImage(mapping[im], voxelsize=im.voxelsize)
            
        segmented_fused_file=timeNamed(segmented_fused_files,t)
        imsave(segmented_fused_file, im)

    
def get_volumes(n, vol, lin_tree):
    """
    Return the volumes of cell n and its progeny up to division
    n : starting cell
    vol : dictionary of volumes
    lin_tree : lineage tree
    """
    tmp=[vol[n]]
    while len(lin_tree.get(n, ''))==1:
        tmp.append(vol[n])
        n=lin_tree[n][0]
    return tmp

def soon_to_divide(c, tree, reverse_tree, len_min,time_begin,time_end):
    """
    Return if a division is too soon (mother divided to late or cell divided too early)
    c : daughter cell to check
    tree : lineage tree
    reverse_tree : reversed lineage tree
    len_min : minimum length allowed
    time_begin: Starting time point
    time_end : final time of the movie
    """
    cell=deepcopy(tree[reverse_tree[c]])
    cell.remove(c)
    cell=cell[0]
    out=[cell]
    while len(tree.get(cell, []))==1:
        cell=tree[cell][0]
        out.append(cell) 
    cell=reverse_tree[c]
    out2=[cell]
    while reverse_tree.has_key(cell):
        cell=reverse_tree[cell]
        out2.append(cell)
    return (len(out)<len_min and out[-1]/10**4!=time_end) or (len(out2)<len_min and out2[-1]/10**4!=time_begin)





def branches_to_delete(tree, volumes, threshold, reverse_tree, soon=False,ShortLifespan=25,PearsonThreshold=0.9,time_begin=1,time_end=192):
    """
    Return the list of cells that have to be removed sorted from the yougest to the oldest
    tree : lineage tree
    volumes : dictionary of volumes information
    reverse_tree : reversed lineage tree
    soon : True if the cell life span has to be taken into account
    """
    from scipy.stats.stats import pearsonr
    nodes=list(set(tree.keys()).union(set([v for values in tree.values() for v in values])))
    leaves=set(nodes)-set(tree.keys())
    last_time=max(leaves)/10**4
    leaves_to_delete=[l for l in leaves if (((l/10**4)<last_time) or (volumes[l]<threshold))]
    branche_max_val=0
    branche_max=[]
    for l in leaves_to_delete:
        b=[]
        b.append(l)
        while len(tree.get(reverse_tree.get(l, ''), ''))==1:
            l=reverse_tree[l]
            b.append(l)
        b.reverse()
        if (min(b)>branche_max_val):
            if len(b)<ShortLifespan: #Time steps length 
                branche_max=list(b)
                branche_max_val=min(b)
            elif reverse_tree.get(min(b), '')!='':
                n1, n2=tree[reverse_tree[min(b)]]
                tmp1, tmp2=get_volumes(n1, volumes, tree), get_volumes(n2, volumes, tree)
                common_len=min(len(tmp1), len(tmp2))
                P=pearsonr(tmp1[:common_len-1], tmp2[:common_len-1])
                if P[0]<-PearsonThreshold:# and deriv>4*10**4:
                    branche_max=list(b)
                    branche_max_val=min(b)
                elif soon and soon_to_divide(min(b), tree, reverse_tree,ShortLifespan,time_begin,time_end):
                    branche_max=list(b)
                    branche_max_val=min(b)

    branche_max.reverse()
    return branche_max

def get_sisters(lin_tree_out, mother, init, end):
    """
    Build the list of all the progeny of the sister cells of init
    lin_tree_out : lineage tree
    mother : common mother of sisters and init cell
    init : cell id
    end : final time point
    """
    sisters=[s for s in lin_tree_out[mother] if s!=init]
    out=[]
    tmp=sisters
    while tmp!=[] and sisters[0]/10**4<=end:
        out.append(sisters)
        tmp=[]
        for s in sisters:
            if lin_tree_out.has_key(s):
                tmp.extend(lin_tree_out[s])
        sisters=tmp
    return out




def remove_too_little_branches(lin_tree, volumes, threshold=2000, soon=False,ShortLifespan=25,PearsonThreshold=0.9,time_begin=1,time_end=192):
    """
    Build a corrected lineage tree (based on volume correlation and cell life span).
    Return a list of cells to fuse, new volumes and new lineage tree
    lin_tree : lineage tree
    volumes : dictionary of volumes
    threshold : volume low threshold for final cells (in voxels)
    soon : True if the cell life span has to be taken into account
    """
    print 'Process remove too little branches'
    reverse_tree={v:k for k, values in lin_tree.iteritems() for v in values}
    branches=branches_to_delete(lin_tree, volumes, threshold, reverse_tree, soon=soon,ShortLifespan=ShortLifespan,PearsonThreshold=PearsonThreshold,time_begin=time_begin,time_end=time_end)
    lin_tree_out=deepcopy(lin_tree)
    to_fuse={}
    new_volumes=deepcopy(volumes)
    while list(branches)!=[]:
        reverse_tree={v:k for k, values in lin_tree_out.iteritems() for v in values}
        b=branches
        last=b[-1]
        mother=reverse_tree.get(last, '')
        if mother!='':
            tmp=get_sisters(lin_tree_out, mother, last, b[0]/10**4)
            for t in tmp:
                time=t[0]/10**4
                cell=[c for c in b if c/10**4==time]
                to_fuse.setdefault(time, []).append((cell[0], t))
                if len(t)==1 and new_volumes.has_key(t[0]) and new_volumes.has_key(cell[0]):
                    new_volumes[t[0]]=new_volumes[cell[0]]+new_volumes[t[0]]
            i=1
            while (i<=len(b)) and (b[-i]/10**4 in [time_cell[0]/10**4 for time_cell in tmp]):
                n=b[-i]
                sis=tmp[np.argwhere(np.array([time_cell[0]/10**4 for time_cell in tmp])==n/10**4)[0, 0]][0]
                lin_tree_out[reverse_tree[n]].remove(n)
                lin_tree_out.pop(n, None)
                if i+1<=len(b):
                    lin_tree_out.setdefault(sis, []).append(b[-(i+1)])
                    reverse_tree[b[-(i+1)]]=sis
                i+=1
        else:
            for n in b:
                lin_tree_out.pop(n, None)
                reverse_tree.pop(n, None)
        branches=branches_to_delete(lin_tree_out, new_volumes, threshold, reverse_tree, soon=soon,ShortLifespan=ShortLifespan,PearsonThreshold=PearsonThreshold,time_begin=time_begin,time_end=time_end)

    #### FUSE TOO EARLY DIVISIONS
    from scipy.stats.stats import pearsonr
    cells_list=[(i[0],i[1]) for i in lin_tree_out.values() if len(i)==2]
    P={}
    deriv={}
    for n1, n2 in cells_list:
        tmp1, tmp2=get_volumes(n1, new_volumes, lin_tree_out), get_volumes(n2, new_volumes, lin_tree_out)
        common_len=min(len(tmp1), len(tmp2))
        if len(tmp1)>8 and len(tmp2)>8: 
            P[(n1, n2)]=pearsonr(tmp1[:common_len-1], tmp2[:common_len-1])
    to_check=[k for k, v in P.iteritems() if v[0]<-.80] #.80 WHAT ? 
    scores_window={}
    for c1, c2 in to_check:
        tmp1, tmp2=get_volumes(c1, volumes, lin_tree_out), get_volumes(c2, volumes, lin_tree_out)
        scores=[]
        common_len=min(len(tmp1), len(tmp2))
        for i in range(0, common_len-4):
            scores.append(pearsonr(tmp1[i:i+5], tmp2[i:i+5])[0])
        scores_window[(c1, c2)]=np.array(scores)

    been_fused={}
    for c, scores in scores_window.iteritems():
        if (np.array(scores)<-.8).all(): #.8 WHAT ? 
            first_size=len(scores)-1
        else:
            out=scores[1:]-scores[:-1]
            sizes=np.argsort(out)[::-1]
            sizes=np.array([s for s in sizes if scores[s]<-.8])
            first_size=sizes[sizes<len(scores)/2]
            if first_size!=[]:
                first_size=first_size[0]
            else:
                first_size=-1
        c1, c2=c
        for i in range(first_size+1):
            next_c1=lin_tree_out.get(c1, [0])[0]
            next_c2=lin_tree_out.get(c2, [0])[0]
            lin_tree_out.setdefault(c1, []).extend(lin_tree_out.get(c2, []))
            for tmp_c in lin_tree_out[c1]:
                reverse_tree[tmp_c]=c1
            lin_tree_out[reverse_tree[c1]]=[c1]
            lin_tree_out.pop(c2, None)
            time=c1/10**4
            to_fuse.setdefault(time, []).append((c2, [c1]))
            new_volumes[c1]=new_volumes[c1]+new_volumes[c2]
            c1=next_c1
            c2=next_c2
            been_fused[c1]=1
    return lin_tree_out, new_volumes, to_fuse, been_fused