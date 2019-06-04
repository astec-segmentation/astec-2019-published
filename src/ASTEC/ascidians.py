from scipy import ndimage as nd
import numpy as np
from copy import copy
from ImageHandling import imread, imsave
from lineage import *


def dist(x,y):
    """
    Return the distance between two points.
    """
    d = np.array(x)-np.array(y)
    return np.sqrt(np.dot(d,d))

def find_neighbors(c1, c2, ref, image):
    """
    Sort cells c1 and c2 in order of their distance to ref in the image 'image'
    c1 : cell label + time of the image * 10**4
    c2 : cell label + time of the image * 10**4
    ref : cell label + time of the image * 10**4
    image : image containing c1, c2 and ref
    """
    from scipy import ndimage as nd
    COM_c1  = np.mean(np.where(image==c1%10**4), axis=1)
    COM_c2  = np.mean(np.where(image==c2%10**4), axis=1)
    COM_ref = np.mean(np.where(image==ref%10**4), axis=1)
    if dist(COM_c1, ref)<dist(COM_c2, ref):
        return [c1, c2], np.mean(np.where((image==(c1%10**4)) | (image==(ref%10**4))), axis=1)
    else:
        return [c2, c1], np.mean(np.where((image==(c2%10**4)) | (image==(ref%10**4))), axis=1)

def order_2_cells_dividing(c1, c2, image):
    """
    Sort 4 cells (given their mother c1, c2):
            1 : the daughter cell of c1 that is the closest to a daughter cell in c2
            2 : the daughter cell of c2 that is the closest to a daughter cell in c1
            3 : the other daughter of c1
            4 : the other daughter of c2
    c1 : list of cells
    c2 : list of cells
    image : segmentated image at time t(c1) + dt
    """
    from scipy import ndimage as nd
    COM_c1=[]
    COM_c2=[]
    COM_c1.append(np.mean(np.where(image==c1[0]%10**4), axis=1))
    COM_c2.append(np.mean(np.where(image==c2[0]%10**4), axis=1))
    COM_c1.append(np.mean(np.where(image==c1[1]%10**4), axis=1))
    COM_c2.append(np.mean(np.where(image==c2[1]%10**4), axis=1))
    distances=[dist(COM_c1[0], COM_c2[0]), dist(COM_c1[0], COM_c2[1]), dist(COM_c1[1], COM_c2[0]), dist(COM_c1[1], COM_c2[1])]
    configuration=np.argmin(distances)
    if configuration==0:
        return [c1[0], c2[0], c1[1], c2[1]], np.mean(np.where((image==(c1[0]%10**4)) | (image==(c2[0]%10**4))), axis=1)
    elif configuration==1:
        return [c1[0], c2[1], c1[1], c2[0]], np.mean(np.where((image==(c1[0]%10**4)) | (image==(c2[1]%10**4))), axis=1)
    elif configuration==2:
        return [c1[1], c2[0], c1[1], c2[0]], np.mean(np.where((image==(c1[1]%10**4)) | (image==(c2[0]%10**4))), axis=1)
    else:
        return [c1[1], c2[1], c1[0], c2[0]], np.mean(np.where((image==(c1[1]%10**4)) | (image==(c2[1]%10**4))), axis=1)

def sort_to_bary(s, bary, image):
    """
    Sort the daughters cells of s from the closest to the point bary to the farest
    s : list of cells
    bary : 3D point
    image : segmentated image at time t(s) + dt
    """
    COM_c1 = np.mean(np.where(image==s[0]), axis=1)
    COM_c2 = np.mean(np.where(image==s[1]), axis=1)
    d1, d2=dist(COM_c1, bary), dist(COM_c2, bary)
    if d1 < d2:
        return s, d2-d1
    else:
        return [s[1], s[0]], d1-d2

def CalculAnimalPole(A,B,segmentation_corrected_files,naming,lin_tree):
    next=[A, B]
    #print next
    im=imread(timeNamed(segmentation_corrected_files,str(A/10**4)))
    time_barycenter={}
    time_barycenter[A/10**4]=np.mean(np.where((im==A%10**4) | (im==B%10**4)), axis=1)
    del im
    treated=set(next)
    # barycenter naming
    while next!=[]:
        prev_A=naming[next[0]]
        prev_B=naming[next[1]]
        nA=lin_tree.get(next[0], [])
        nB=lin_tree.get(next[1], [])
        print ' Process from '+str(prev_A)+ ' and ' + str(prev_B)
        treated=treated.union(set(next))
        if len(nA)==1 and len(nB)==1:
            naming[nA[0]] = naming[next[0]]
            naming[nB[0]] = naming[next[1]]
            im=imread(timeNamed(segmentation_corrected_files,str(nA[0]/10**4)))
            time_barycenter[nA[0]/10**4]=np.mean(np.where((im==nA[0]%10**4) | (im==nB[0]%10**4)), axis=1)
            print 'at '+str(nA[0]/10**4)+' ->' +str(time_barycenter[nA[0]/10**4])
            del im
            next=[nA[0], nB[0]]
        elif len(nA)==1 and len(nB)==2:
            #Recupere la plus proche du barycentre
            naming[nA[0]]=naming[next[0]]
            nB, time_barycenter[nA[0]/10**4]=find_neighbors(nB[0], nB[1], nA[0],
                imread(timeNamed(segmentation_corrected_files,str(nA[0]/10**4))))
            for i, c in enumerate(nB):
                naming[c] = 'b' + str(int(prev_B.split('.')[0][1:])+1) + '.' + \
                    '%04d'%(int(prev_B.split('.')[1][:-1])*2-1+i) + '*'
            print 'at '+str(nA[0]/10**4)+' ->' +str(time_barycenter[nA[0]/10**4])
            next=[nA[0], nB[0]]

        elif len(nB)==1 and len(nA)==2:
            naming[nB[0]]=naming[next[1]]
            nA, time_barycenter[nA[0]/10**4]=find_neighbors(nA[0], nA[1], nB[0],
                imread(timeNamed(segmentation_corrected_files,str(nA[0]/10**4))))
            for i, c in enumerate(nA):
                naming[c] = 'a' + str(int(prev_A.split('.')[0][1:])+1) + '.' + \
                    '%04d'%(int(prev_A.split('.')[1][:-1])*2-1+i) + '_'
            print 'at '+str(nA[0]/10**4)+' ->' +str(time_barycenter[nA[0]/10**4])
            next=[nA[0], nB[0]]

        elif len(nB)==len(nA)==2:
            next=[]
            order, time_barycenter[nA[0]/10**4]=order_2_cells_dividing(nA, nB,
                imread(timeNamed(segmentation_corrected_files,str(nA[0]/10**4))))
            for i, c in enumerate(order):
                if i%2:
                    prev=prev_B
                else:
                    prev=prev_A
                naming[c] = prev[0] + str(int(prev.split('.')[0][1:])+1) + '.' + \
                    '%04d'%(int(prev.split('.')[1][:-1])*2-1+i) + prev[-1]
            print 'at '+str(nA[0]/10**4)+' ->' +str(time_barycenter[nA[0]/10**4])
            next=[order[0], order[1]]

        else:
            print "Finish "
            next=[]


    return time_barycenter,treated

def name_propagation(lin_tree, starting_name, segmentation_corrected_files, begin,end, step):
    """
    Return the propagated names of each cell in the lineage tree using Conklin rule (eulerian distand, not geodesic)
    lin_tree : lineage tree
    starting_name : manually processed starting names ({cell id : name}) the names has to be this format : 'a7.xxxx*' or 'a7.xxxx_'
    path : path toward the segmentated image (the names of the segmented images has to be fw_seg_t002_corrected.inr, TO CHANGE)
    end : ending time point
    step : dt value
    """
    print 'Process Name Propagation'
    from copy import copy
    inv_lin_tree={ d : m for m, daughters in lin_tree.iteritems() for d in daughters }
    naming={}
    for k, v in starting_name.iteritems():
        naming[k]=v
 

    stage=naming.values()[0].split('.')[0][1:]
    #First we calcul the Animale Pole to orient named cell division 
    A=naming.keys()[naming.values().index("a"+stage+".0001*")]
    B=naming.keys()[naming.values().index("b"+stage+".0001_")]
    time_barycenter1,treated1=CalculAnimalPole(A,B,segmentation_corrected_files,naming,lin_tree)
    A=naming.keys()[naming.values().index("a"+stage+".0001_")]
    B=naming.keys()[naming.values().index("b"+stage+".0001*")]
    time_barycenter2,treated2=CalculAnimalPole(A,B,segmentation_corrected_files,naming,lin_tree)
    treated=treated1.union(treated2)
    #Add starting name
    for cell in naming:
        if cell not in treated:
            treated.add(cell)

    #Average Barycenters
    time_barycenter={}
    for t in range(begin, end+1, step):
        n=0
        coord=np.zeros(3)
        if t in time_barycenter1:
            coord+=time_barycenter1[t]
            n+=1
        if t in time_barycenter2:
            coord+=time_barycenter2[t]
            n+=1
        if n==0:
            time_barycenter[t]=time_barycenter[t-1]
        else:
            time_barycenter[t]=np.divide(coord,n)
        #print ' at '+str(t)+ '->'+str(time_barycenter[t])


    error_metric={}


    pre_treated=np.array(list(treated))
    for i in range(begin, end+1, step):
        print "Progagate cell name at "+str(i)
        im = imread(timeNamed(segmentation_corrected_files,i))
        ids=i*10**4+np.unique(im[im>1]).astype(np.uint32)
        treated=list(pre_treated[pre_treated/10**4==i])
        for cell in ids:
            print " Process cell "+str(cell)
            if not cell in treated:
            	#print ' -> not treated'
                treated.append(cell)
                mother=inv_lin_tree.get(cell, None)
                sibling=lin_tree.get(mother, cell)
                if sibling==[cell]:
                    naming[cell]=naming.get(mother, "")
                elif cell in sibling and len(sibling)==2:
                    #print ' Sinling sibling[0]%10**4='+str(sibling[0]%10**4) 
                    #print ' Sinling sibling[1]%10**4='+str(sibling[1]%10**4) 
                    #print time_barycenter
                    #print i
                    #print time_barycenter[i]
                    cells_sorted, dist_diff=sort_to_bary([sibling[0]%10**4, sibling[1]%10**4], time_barycenter[i], im)
                    error_metric[cell]=dist_diff
                    for j, c in enumerate(cells_sorted):
                        prev_n=naming.get(mother, "")
                        if prev_n:
                            naming[c+i*10**4] = prev_n[0] + str(int(prev_n.split('.')[0][1:])+1) + '.' + \
                    '%04d'%(int(prev_n.split('.')[1][:-1])*2-1+j) + prev_n[-1]
                        treated.append(c+i*10**4)
                else:
                    treated.extend(sibling)
            if not cell in naming:
                naming[cell]=""
            print " -->" +naming[cell]
  
    """
    #suite a erreurs sur le naming de Ralph
    #Propage Errors
    errors=error_metric
    roots=[k for k in lin_tree.keys() if k<=2*10**4]
    todo=roots
    inv_tree={ v : k for k, values in lin_tree.iteritems() for v in values }
    while todo!=[]:
        n=todo[0]
        if len(lin_tree.get(n, ''))==2:
            tmp=lin_tree[n]
            if errors.get(tmp[0], 0)==0:
                errors[tmp[0]]=errors[tmp[1]]
            else:
                errors[tmp[1]]=errors[tmp[0]]
        if errors.get(n, '')=='':
            errors[n]=errors.get(inv_tree.get(n, ''), 20)
        todo.remove(n)
        new=lin_tree.get(n, [])
        todo.extend(new)

    return naming,errors
    """
    return naming



def find_color(name, cell_fate):
    for k, v in cell_fate.iteritems():
        if name[:-1].replace(' ', '') in v[0]:
            return k
    return ''

def fate_map(lin_tree,names,begin,end,step):
    print 'Process Fate Map'

    cell_fate2={
        'Anterior Endoderm': (["a7.0001", "a7.0002", "a7.0005"], 1),
        'Posterior Endoderm':(["b7.0001", "b7.0002", "b9.0034"], 2),

        'germ line': (["b7.0006"], 3),

        'Mesoderm 1 Notochord':(["a7.0003", 'a7.0007'], 4),
        'Mesoderm 2 Notochord':(["b8.0006"], 5),
        'Mesoderm Trunk Lateral Cell':(['a7.0006'], 6),
        'Mesoderm Trunk ventral Cell':(['b7.0005'], 7),
        'Mesoderm First Muscle':(['b7.0004', 'b7.0008'], 8),
        'Mesoderm Second Muscle':(['a9.0031', 'b9.0033'], 9),
        'Mesoderm Mesenchyme':(['b7.0007', 'b8.0005'], 10),

        'Posterior ventral Neural Plate': (["a7.0004"],11),
        'Anterior + Dorsal Neural Plate': (['a7.0009', 'a7.0010', 'b8.0019'], 12),
        'Lateral Neural Plate':(['a7.0013', 'a8.0015', 'a9.0032'], 13),
        
        'Trunk Epidermis':(['a7.0011', 'a7.0012', 'a7.0014', 'a7.0015', 'a7.0016'], 14),
        'Midline Tail Epidermis':(['b8.0020', 'b8.0018', 'b9.0041', 'b8.0027', 'b9.0056', 'b9.0062', 'b9.0064'], 15),
        'Mediolateral Tail Epidermis':(['b8.0024', 'b9.0045', 'b9.0042', 'b9.0043', 'b9.0049', 'b9.0055', 'b9.0061', 'b9.0063'], 16),
        'Lateral Tail Epidermis':(['b9.0044', 'b8.0026', 'b9.0050', 'b9.0046', 'b8.0029', 'b8.0030'], 17)
    }

    cell_fate={
        'Endoderm': (["a7.0001", "a7.0002", "a7.0005", "b7.0001", "b7.0002", "b9.0034"], 7),
        'germ line': (["b7.0006"], 3),
        'Epidermis': (['a7.0011', 'a7.0012', 'a7.0014', 'a7.0015', 'a7.0016', 'b8.0020', 'b8.0018', 'b7.0011', 'b7.0012', 'b7.0013', 'b7.0014', 'b7.0015', 'b7.0016'], 2),
        'Mesoderm': (["a7.0003", 'a7.0006','a7.0007', "b7.0003", 'b7.0005', 'b7.0004', 'b7.0008', 'a9.0031', 'b9.0033', 'b7.0007'], 6),
        'Nervous system': (["a7.0004", "a7.0009", "a7.0010", 'a7.0013', 'a8.0015', 'a9.0032', 'b8.0019'], 5)
    }

    inv_lin_tree={ v : k for k, values in lin_tree.iteritems() for v in values }
    fate={}
    nodes=list(set(lin_tree.keys()).union(set([v for values in lin_tree.values() for v in values])))
    name_fate={}
    for t in range(begin, end+1, step):
        for n in nodes:
            if t*10**4<n<(t+1)*10**4: 
                if fate.get(inv_lin_tree.get(n, ''), '')!='':
                    fate[n]=fate[inv_lin_tree[n]]
                    name_fate[names.get(n, '')]=fate[n]
                else:
                    col=find_color(names.get(n, ''), cell_fate)
                    if col!='':
                        fate[n]=col
                        name_fate[names.get(n, '')]=col

    fate2={}
    for t in range(begin, end+1, step):
        for n in nodes:
            if t*10**4<n<(t+1)*10**4: 
                if fate2.get(inv_lin_tree.get(n, ''), '')!='':
                    fate2[n]=fate2[inv_lin_tree[n]]
                else:
                    col=find_color(names.get(n, ''), cell_fate2)
                    if col!='':
                        fate2[n]=col

    return fate,fate2
   


#Surface of Contacts
def __slices_dilation(slices, maximum=[np.inf, np.inf, np.inf]):
    return tuple([slice(max(0, s.start-1), min(s.stop+1, maximum[i])) for i, s in enumerate(slices)])

def slices_dilation(slices, maximum=[np.inf, np.inf, np.inf], iterations=1):
    for i in range(iterations):
        slices=__slices_dilation(slices, maximum)
    return slices    

def process_contact_surface(t, seg_file_format):
    seg=imread(timeNamed(seg_file_format,t))
    bboxes=nd.find_objects(seg)
    cells=np.unique(seg)
    contact_surf={}
    seg[seg==0]=1
    labels=cells[2:]
    barycenters_tmp=dict(zip(t*10**4+labels, nd.center_of_mass(np.ones_like(seg), seg, labels)))
    barycenters={k:(v[0]*2, v[1]*2, v[2]*2) for k, v in barycenters_tmp.iteritems()}
    for c in cells:
        if c!=1 and c!=0:
            bb=slices_dilation(bboxes[c-1], iterations=2)
            tmp=seg[bb]
            bounds=tmp[(nd.binary_dilation(tmp==c) & (tmp!=c))]
            contact_surf[c]=dict(zip(np.unique(bounds), nd.sum(np.ones_like(bounds), bounds, np.unique(bounds))))
    surfaces={}
    for k, v in contact_surf.iteritems():
        surfaces[t*10**4+k]={t*10**4+key:val for key, val in v.iteritems()}
    return contact_surf, surfaces, barycenters
        

def process_multiple_contacts(seg_file_format, begin, end, delta):
    print 'Process Surface of Contact'
    contact_surf, surfaces, barycenters, inertia_axis, eig_values = {}, {}, {}, {}, {}
    for t in range(begin, end+1, delta):
        print ' ->Surface of contact at '+str(t)
        cs, s, b = process_contact_surface(t, seg_file_format)
        contact_surf.update(cs)
        surfaces.update(s)
        barycenters.update(b)

    return surfaces,barycenters,contact_surf
