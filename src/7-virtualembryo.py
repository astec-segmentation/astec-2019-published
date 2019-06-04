#!/usr/bin/python2.7



#Fusion Process
from definitions import *
if not os.path.isdir(vo_Path):
    os.mkdir(vo_Path)  
os.system("cp -f "+astec_Path+"definitions.py "+vo_Path )
os.system("cp -f "+astec_Path+"7-virtualembryo.py "+vo_Path )

from ImageHandling import imread, imsave, SpatialImage
from lineage import read_lineage_tree,timeNamed

sys.path.append(vo_lib_Path)
from oascripts import image_to_vtk_cell,vtk_polydata_to_cell_triangular_meshes,composed_triangular_mesh,save_obj_to_cells,loadOBJ



#CREATE MESH OBJ FILE

#PARAMETERS
mesh_fineness=0.5;
decimate="02"
decimateLow="2"
sub_factor=1
decimator=vo_lib_Path+'MeshLabResamble' + str(decimate) + '.mlx'
decimatorLow=vo_lib_Path+'MeshLabResamble' + str(decimateLow) + '.mlx'
smother=vo_lib_Path+'Taubin0.5--0.3-10.mlx';
meshlabserver='meshlabserver';

os.system('Xvfb :100 &')
display="export DISPLAY=:100.0 " #DONT FORGET TO ADD THAT


for t in range(begin,end+1): #END
    #segmented_img=imread('TEST.tif')#segmented_img = img = imread(timeNamed(postsegment_files,t))
    vo_file=timeNamed(vo_files,t)
    if not os.path.isfile(vo_file):
        segmented_img = img = imread(timeNamed(postsegment_files,t))
        cell_path=vo_Path+'TEMP_'+str(t)+'/'
        if not os.path.isdir(cell_path):
            os.mkdir(cell_path)  
        print 'Generate Cells as OBJ in  ' +cell_path
        cell_polydata = image_to_vtk_cell(segmented_img,coef=1.0,mesh_fineness=mesh_fineness)
        cell_triangular_meshes = vtk_polydata_to_cell_triangular_meshes(cell_polydata)
        mesh, matching = composed_triangular_mesh(cell_triangular_meshes)
        save_obj_to_cells(mesh, matching,cell_path)
        #DECIMATE AND MERGE CELLS IN OBJ
        vo_file=timeNamed(vo_files,t)
        fLQ = open(vo_file, "w") 
        shiftnLQ=0
        for filename in os.listdir(cell_path):
            if filename.find("cell")==0 and  filename.find("smooth")==-1 and  filename.find("decimate")==-1:
                cellid=int(filename.replace("cell_","").replace(".obj",""))
                print 'Smooth cell ' +str(cellid)
                smoothname=cell_path+filename.replace(".obj","_smooth.obj")
                os.system(meshlabserver +' -i '+cell_path+filename+' -o '+smoothname+' -s ' + smother);
                print 'Decimate cell ' +str(cellid)
                decimatefilename=cell_path+filename.replace(".obj","_decimate.obj")
                os.system(meshlabserver +' -i '+smoothname+' -o '+decimatefilename+' -s '+ decimator+' ');
                if not os.path.isfile(decimatefilename):
                	os.system(meshlabserver +' -i '+smoothname+' -o '+decimatefilename+' -s '+ decimatorLow+' ');
                vertex,face=loadOBJ(decimatefilename)
                fLQ.write('g cell_' +str(cellid) +"\n")
                for v in vertex:
                    fLQ.write('v ' + str(v[0]) +' '+str(v[1]) +' '+str(v[2])+'\n')
                    #fLQ.write('v ' + str(v[0]-CenterEmbryo[0]) +' '+str(v[1]-CenterEmbryo[1]) +' '+str(v[2]-CenterEmbryo[2])+'\n')
                for fs in face:
                    fLQ.write('f ' + str(fs[0]+shiftnLQ) +' '+str(fs[1]+shiftnLQ) +' '+str(fs[2]+shiftnLQ)+'\n')
                shiftnLQ=shiftnLQ+len(vertex)
        fLQ.close()
        os.system('rm -rf '+cell_path)




#UPLOAD IT IN THE DB VIRTUAL EMBRYO
from VirtualEmbryo import VirtualEmbryo
ve=VirtualEmbryo(vo_login,vo_passwd)



#CREATE DATA
#il faut tester d'abord si le dataset n'existe pas
dataset=ve.selectDataSetByName(EN)
if ve.isDataSet():
    print ' Data set '+EN+' already exist'
    if not ve.ownDataSet():
        print " but you cannot change anything on this embryo, you're not the owner..."
        quit()
    else:
        print ' -> delete previous dataset before upload'
        ve.deleteDataSet()
ve.createDataSet(EN,begin,end)

#UPLOAD OBJ
for t in range(begin,end+1):
    vo_file=timeNamed(vo_files,t)
    print 'Upload cells at '+str(t)
    f=open(vo_file,'r')
    #We have to convert first face index
    idx=0
    newidx=0
    obj=""
    for line in f:
        if line.find('g ')==0: # NEW CELL
            idx=newidx
            cell=line.split(" ")[1]
            obj+=line
        if line.find('v ')==0: # VERTEX
            obj+=line
        if line.find('f ')==0: # FACE
            faces=line.split(" ")
            obj+='f'
            for i in range(1,4):
                obj+=' '+str(int(faces[i])-idx)
                newidx=max(newidx,int(faces[i]))
            obj+='\n'
    f.close()
    ve.timepoint[t]=obj
    ve.uploadTimePoint(t)

#Read Lineage Infos
lin_tree_information=read_lineage_tree(named_lineage_tree_filename)


def getID( longID):
    t = round(longID / 10000);
    return int(longID - t * 10000);
def getTime ( longID):
    return  int(round(longID / 10000));

#Upload Lineage
lineage=lin_tree_information['lin_tree']
lineageText=""
for mother, daughters in lineage.iteritems():
    for daughter in daughters:
        lineageText+=str(getTime(mother)) +';'+str(getID(mother)) +';'+ str(getID(daughter))+'\n'
ve.uploadInfos("Lineage",lineageText,'string')

#Upload Volumes
volume=lin_tree_information['volumes_information']
volumeText=""
for key, value in volume.iteritems():
    volumeText+=str(key) +':'+str(value)+'\n'
ve.uploadInfos("Volume",volumeText,'float')

#Upload H Mins
hmins=lin_tree_information['h_mins_information']
hminsText=""
for key, value in hmins.iteritems():
    hminsText+=str(key) +':'+str(value)+'\n'
ve.uploadInfos("H Min",hminsText,'float')

#Upload Sigmas
sigmas=lin_tree_information['sigmas_information']
sigmasText=""
for key, value in sigmas.iteritems():
    sigmasText+=str(key) +':'+str(value)+'\n'
ve.uploadInfos("Sigmas",sigmasText,'float')


#Upload Names
Names=lin_tree_information['Names']
NamesText=""
for key, value in Names.iteritems():
    NamesText+=str(key) +':'+str(value)+'\n'
ve.uploadInfos("Name",NamesText,'string')


#Upload Neigbhors with surface of contact
ccsurface=lin_tree_information['surface']
ccsurfaceText=""
for key, value in ccsurface.iteritems():
    ccel=str(key) +':{'
    for keyX,valueX in value.iteritems():
        ccel+=str(keyX) +':'+str(valueX)+','
    if ccel[-1]==',':
        ccel=ccel[:len(ccel)-1]
    ccel+='}\n'
    ccsurfaceText+=ccel
ve.uploadInfos("Neighbors",ccsurfaceText,'dict')


#Upload Fate
fate=lin_tree_information['fate']
fateText=""
for key, value in fate.iteritems():
    fateText+=str(key) +':'+str(value)+'\n'
ve.uploadInfos("Fate",fateText,'string')

#Create a fatemap selection
Fates=[]
for key, value in fate.iteritems():
    if value.strip() not in Fates:
        Fates.append(value.strip())

print Fates
def getFateIdx(Fates,fate):
    i=1
    for f in Fates:
        if f==fate:
            return i
        i+=1
    return -1

ve.createSelection("FateMap")
for key, value in fate.iteritems():
    ve.addSelectionToCell(key,getFateIdx(Fates,value.strip()))
ve.updateSelection() 

