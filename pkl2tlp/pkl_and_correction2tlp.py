#!/usr/bin/env/python

###########################################################################
###########################################################################
## Copyright (C) 2018  Guignard Leo <guingardl__@__janelia.hhmi.org>     ##
##                                                                       ##
## This program is free software: you can redistribute it and/or modify  ##
## it under the terms of the GNU General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or     ##
## (at your option) any later version.                                   ##
##                                                                       ##
## This program is distributed in the hope that it will be useful,       ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of        ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         ##
## GNU General Public License for more details.                          ##
##                                                                       ##
## You should have received a copy of the GNU General Public License     ##
## along with this program.  If not, see <http://www.gnu.org/licenses/>. ##
###########################################################################
###########################################################################

import numpy as np
from copy import deepcopy
import argparse
import cPickle as pkl
from matplotlib import pyplot as plt

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




def remove_too_little_branches(lin_tree, volumes, threshold=2000, soon=False,ShortLifespan=25,PearsonThreshold=0.9,
                                PTh2=.8,time_begin=1,time_end=192):
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
    branches=branches_to_delete(lin_tree, volumes, threshold, reverse_tree, soon=soon,ShortLifespan=ShortLifespan,
                                PearsonThreshold=PearsonThreshold,time_begin=time_begin,time_end=time_end)
    lin_tree_out=deepcopy(lin_tree)
    to_fuse={}
    new_volumes=deepcopy(volumes)

    been_fused={}
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
                been_fused[cell[0]] = 2
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
        branches=branches_to_delete(lin_tree_out, new_volumes, threshold, reverse_tree, soon=soon,
                                    ShortLifespan=ShortLifespan,PearsonThreshold=PearsonThreshold,time_begin=time_begin,time_end=time_end)

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
    to_check=[k for k, v in P.iteritems() if v[0]<-PTh2] #.80 WHAT ? 
    scores_window={}
    for c1, c2 in to_check:
        tmp1, tmp2=get_volumes(c1, volumes, lin_tree_out), get_volumes(c2, volumes, lin_tree_out)
        scores=[]
        common_len=min(len(tmp1), len(tmp2))
        for i in range(0, common_len-4):
            scores.append(pearsonr(tmp1[i:i+5], tmp2[i:i+5])[0])
        scores_window[(c1, c2)]=np.array(scores)

    for c, scores in scores_window.iteritems():
        if (np.array(scores)<-PTh2).all(): #.8 WHAT ? 
            first_size=len(scores)-1
        else:
            out=scores[1:]-scores[:-1]
            sizes=np.argsort(out)[::-1]
            sizes=np.array([s for s in sizes if scores[s]<-PTh2])
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

def write_tlp_from_lin_tree(name, lin_tree_information, ListProperties = None, extra_property = None):
    """
    Write a lineage tree into an understable tulip file
    name : path to the tulip file to create
    lin_tree : lineage tree to write
    properties : dictionary of properties { 'Property name': [{c_id: prop_val}, default_val]}
    """
    
    lin_tree=lin_tree_information['lin_tree']
    if ListProperties is None:
        ListProperties=['volumes_information','h_mins_information','sigmas_information']
    properties=[] #Properties Dictionary initialisation
    if ListProperties is None:
        ListProperties = []
    if extra_property is None:
        extra_property = []
    for l in ListProperties:
        if l is not 'Names':
            properties.append((l,lin_tree_information[l],sum(lin_tree_information[l].values())/float(len(lin_tree_information[l]))))


    nodes=set(lin_tree.keys()).union(set([v for values in lin_tree.values() for v in values]))

    f=open(name, "w")
    inv_lin_tree={v:k for k, vals in lin_tree.iteritems() for v in vals}
    f.write("(tlp \"2.0\"\n")
    f.write("(nodes ")
    for n in nodes:
        f.write(str(n)+ " ")
    f.write(")\n")


    nodes=set(lin_tree.keys()).union(set([v for values in lin_tree.values() for v in values]))
    count_edges=0
    for m, ds in lin_tree.iteritems():
        count_edges+=1
        for d in ds:
            f.write("(edge " + str(count_edges) + " " + str(m) + " " + str(d) + ")\n")
    f.write("(property 0 int \"id\"\n")
    f.write("\t(default \"0\" \"0\")\n")
    for node in nodes:
        f.write("\t(node " + str(node) + str(" \"") + str(node) + "\")\n")
    f.write(")\n")

    for property in properties:
        prop_name=property[0]
        vals=property[1]
        default=property[2]
        f.write("(property 0 string \""+prop_name+"\"\n")
        f.write("\t(default \""+str(default)+"\" \"0\")\n")
        for node in nodes:
            f.write("\t(node " + str(node) + str(" \"") + str(vals.get(node, default)) + "\")\n")
        f.write(")\n") 
    for property in extra_property:
        prop_name=property[0]
        vals=property[1]
        default=property[2]
        f.write("(property 0 double \""+prop_name+"\"\n")
        f.write("\t(default \""+str(default)+"\" \"0\")\n")
        for node in nodes:
            f.write("\t(node " + str(node) + str(" \"") + str(vals.get(node, default)) + "\")\n")
        f.write(")\n") 

    f.write(")")
    f.close()


def main():
    parser = argparse.ArgumentParser(description='Convert pkl lineage into tulip lineage.')
    parser.add_argument('-i', '--input', help='input pickle .pkl file', required=True)
    parser.add_argument('-o', '--output', help='output tulip file (has to end with .tlp)', required=True)
    parser.add_argument('-PT1', '--PearsonTh1', help='First Pearson Threshold (for the stopping branches), between 0 and 1, default 0.9',
                                                default=0.9, type=float)
    parser.add_argument('-PT2', '--PearsonTh2', 
                                help='Second Pearson Threshold (for the divisions that happen too early), between 0 and 1, default 0.8', 
                                default=0.8, type=float)
    parser.add_argument('-SLS', '--ShortLifespan', 
                                help='Branch length bellow which cells are considered too short to be true, default 25', 
                                default=25, type=int)

    
    args = parser.parse_args()
    with open(args.input) as f:
        DATA = pkl.load(f)

    lineage = DATA['lin_tree']
    volumes = DATA['volumes_information']

    lin_tree_out, new_volumes, to_fuse, been_fused = remove_too_little_branches(lineage, volumes,
                                                        ShortLifespan=args.ShortLifespan,
                                                        PearsonThreshold=args.PearsonTh1,
                                                        PTh2=args.PearsonTh2)

    write_tlp_from_lin_tree(args.output.replace('.tlp', '_before.tlp'), DATA, extra_property = [['to fuse', been_fused, 0]])


    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1)
    cells_before = set([c for c in lin_tree_out] + [ci for c in lin_tree_out.itervalues() for ci in c])
    times = [c//10**4 for c in cells_before]
    times = np.unique(times)

    cells = set([c for c in lin_tree_out] + [ci for c in DATA['lin_tree'].itervalues() for ci in c])
    cell_per_TP = np.array([0, ]*len(times))
    for c in cells:
        cell_per_TP[np.where(times == c//10**4)[0]] += 1
    ax.plot(times, cell_per_TP, lw = 2, alpha = .7, label = 'Before')

    cell_per_TP = np.array([0, ]*len(times))
    for c in cells_before:
        cell_per_TP[np.where(times == c//10**4)[0]] += 1

    ax.plot(times, cell_per_TP, lw = 2, alpha = .7, label = 'Corrected')

    ax.set_xlabel('Time [TP]', size = 26)
    ax.set_ylabel('#cells', size = 26)
    ax.legend(fontsize = 22, loc = 'upper left')
    ax.tick_params(labelsize = 20)
    fig.tight_layout()
    fig.savefig(args.output.replace('.tlp', '.pdf'))
    plt.close('all')

    DATA['lin_tree'] = lin_tree_out
    DATA['volumes_information'] = new_volumes

    write_tlp_from_lin_tree(args.output, DATA)


if __name__ == '__main__':
    main()