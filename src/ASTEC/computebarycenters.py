#Fusion Process
from definitions import *

from ImageHandling import imread, imsave, SpatialImage
from lineage import read_lineage_tree,timeNamed
import numpy as np 
import cPickle as pkl
from threading import Thread

def save_pkl(obj, filename):
    #saves pickle object, with name given by "filename".
    with open(filename, 'wb') as output:
        pkl.dump(obj, output, pkl.HIGHEST_PROTOCOL)

def load_pkl(filename):
    #loads pickle object
    with open(filename, 'r') as file:
        obj=pkl.load(file)
    return(obj)



class Thread_BaryCenter(Thread):
    def __init__(self,filename,t):
        Thread.__init__(self) 
        self.filename=filename
        self.t=t
 
    def run(self):
        print 'Read '+self.filename
        baryT={}
        segmented_img = imread(self.filename)
        cells=np.unique(segmented_img)
        cells=cells[np.where(cells!=1)] #1 IS BACKGROUND
        for cell in cells:
            print 'Process cell'+str(cell)+' at '+str(self.t)
            [x,y,z]=np.where(segmented_img==cell)
            baryT[self.t*10**4+cell]=[x.mean(),y.mean(),z.mean()]
        save_pkl(baryT,glastec_postsegment_Path+'new_barycenters_t'+str(self.t)+'.pkl')

def WaitThreadsEnd(Threads):
    nb=0
    for th in Threads:
        print ' Wait thread ' + str(nb) + ' / '+str(len(Threads))
        th.join() #Wait the end
        nb+=1
   


